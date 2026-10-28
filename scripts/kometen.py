#!/usr/bin/env python3
"""
Lädt die Bahnelemente aller Kometen vom Minor Planet Center (MPC) und schreibt
comets.json mit den Kometen, die in den nächsten 60 Tagen heller als 13,5 mag
werden könnten. Läuft täglich als GitHub Action, braucht keine Zusatzpakete.
"""
import json, math, time, urllib.request

URL = "https://www.minorplanetcenter.net/iau/MPCORB/CometEls.txt"
MAG_LIMIT = 13.5
K = 0.01720209895          # Gaußsche Gravitationskonstante
D = math.pi / 180
EPS = 23.4392911 * D       # Schiefe der Ekliptik J2000


def jd_from_date(y, m, d):
    if m <= 2:
        y -= 1; m += 12
    a = y // 100; b = 2 - a + a // 4
    return math.floor(365.25 * (y + 4716)) + math.floor(30.6001 * (m + 1)) + d + b - 1524.5


def helio(c, jd):
    """Heliozentrische ekliptikale Koordinaten (J2000) und Sonnenabstand."""
    q, e, dt = c["q"], c["e"], jd - c["T"]
    if e < 0.999999:
        a = q / (1 - e); n = K / a ** 1.5
        M = math.remainder(n * dt, 2 * math.pi)
        E = M + 0.85 * e * (1 if math.sin(M) >= 0 else -1) if e > 0.8 else M
        for _ in range(100):
            dE = (E - e * math.sin(E) - M) / (1 - e * math.cos(E)); E -= dE
            if abs(dE) < 1e-12: break
        r = a * (1 - e * math.cos(E))
        nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(E / 2), math.sqrt(1 - e) * math.cos(E / 2))
    elif e <= 1.000001:
        W = 3 * K / math.sqrt(2 * q ** 3) * dt
        Y = (W / 2 + math.sqrt(W * W / 4 + 1)) ** (1 / 3)
        s = Y - 1 / Y
        nu = 2 * math.atan(s); r = q * (1 + s * s)
    else:
        a = q / (e - 1); N = K * dt / a ** 1.5
        H = math.asinh(N / e)
        for _ in range(100):
            dH = (e * math.sinh(H) - H - N) / (e * math.cosh(H) - 1); H -= dH
            if abs(dH) < 1e-12: break
        r = a * (e * math.cosh(H) - 1)
        nu = 2 * math.atan(math.sqrt((e + 1) / (e - 1)) * math.tanh(H / 2))
    w, O, i = c["w"] * D, c["node"] * D, c["i"] * D
    u = w + nu
    x = r * (math.cos(O) * math.cos(u) - math.sin(O) * math.sin(u) * math.cos(i))
    y = r * (math.sin(O) * math.cos(u) + math.cos(O) * math.sin(u) * math.cos(i))
    z = r * math.sin(u) * math.sin(i)
    return x, y, z, r


def sun_geo(jd):
    n = jd - 2451545.0
    L = (280.460 + 0.9856474 * n) % 360
    g = ((357.528 + 0.9856003 * n) % 360) * D
    lam = (L + 1.915 * math.sin(g) + 0.020 * math.sin(2 * g) - 1.397 * n / 36525) * D
    R = 1.00014 - 0.01671 * math.cos(g) - 0.00014 * math.cos(2 * g)
    return R * math.cos(lam), R * math.sin(lam), 0.0


def observe(c, jd):
    x, y, z, r = helio(c, jd)
    sx, sy, sz = sun_geo(jd)
    X, Y, Z = x + sx, y + sy, z + sz
    delta = math.sqrt(X * X + Y * Y + Z * Z)
    ye = Y * math.cos(EPS) - Z * math.sin(EPS)
    ze = Y * math.sin(EPS) + Z * math.cos(EPS)
    ra = math.degrees(math.atan2(ye, X)) % 360
    dec = math.degrees(math.asin(ze / delta))
    mag = c["H"] + 5 * math.log10(delta) + 2.5 * c["G"] * math.log10(r)
    return ra, dec, mag, r, delta


def num(s):
    s = s.strip()
    return float(s) if s else None


def parse(line):
    """Ein Datensatz im MPC-Exportformat für Kometenbahnen."""
    if len(line) < 105:
        return None
    t = line[14:91].split()
    if len(t) < 8:
        return None
    H, G = num(line[91:95]), num(line[96:100])
    if H is None:
        return None
    y, m, d = int(t[0]), int(t[1]), float(t[2])
    return {
        "name": line[102:158].strip(),
        "T": round(jd_from_date(y, m, d), 5),
        "q": float(t[3]), "e": float(t[4]), "w": float(t[5]),
        "node": float(t[6]), "i": float(t[7]),
        "H": H, "G": G if G is not None else 4.0,
    }


def main():
    req = urllib.request.Request(URL, headers={"User-Agent": "deepsky-planer (GitHub Action)"})
    text = urllib.request.urlopen(req, timeout=60).read().decode("utf-8", "replace")
    now = time.time() / 86400 + 2440587.5
    keep = []
    for line in text.splitlines():
        try:
            c = parse(line)
        except (ValueError, IndexError):
            continue
        if not c or c["q"] > 8:
            continue
        try:
            obs = [observe(c, now + dd) for dd in (0, 10, 20, 30, 45, 60)]
        except (ValueError, ZeroDivisionError, OverflowError):
            continue
        best = min(o[2] for o in obs)
        if best <= MAG_LIMIT and max(o[1] for o in obs) > -45:
            c["maxMag60"] = round(best, 1)
            keep.append(c)
    keep.sort(key=lambda c: c["maxMag60"])
    out = {"updated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "source": "Minor Planet Center, CometEls.txt", "comets": keep}
    with open("comets.json", "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print(f"{len(keep)} Kometen gespeichert")


if __name__ == "__main__":
    main()

