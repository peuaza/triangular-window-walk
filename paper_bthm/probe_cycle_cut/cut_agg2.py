import json

d = json.load(open("cut_cells.json", encoding="utf-8"))
R = d["rules"]
NS = ("512", "1024", "2048")
KEYS = ("tapes", "b0_meet", "b0_cyc", "bC_meet", "bC_rescue", "bC_cap", "bC_empty")


def tot(n):
    s = {k: 0 for k in KEYS}
    for r, by in R.items():
        c = by.get(n)
        if not c:
            continue
        for k in KEYS:
            s[k] += c.get(k) or 0
    return s


for n in NS:
    s = tot(n)
    print(n, "tapes", s["tapes"],
          "meet0 %.4f" % (s["b0_meet"] / s["tapes"]),
          "cap0 %.4f" % (s["b0_cyc"] / s["tapes"]),
          "meetC %.4f" % (s["bC_meet"] / s["tapes"]),
          "rescued %d of %d %.4f" % (s["bC_rescue"], s["b0_cyc"],
                                     (s["bC_rescue"] / s["b0_cyc"]) if s["b0_cyc"] else 0),
          "capC", s["bC_cap"])

T = sum(tot(n)["tapes"] for n in NS)
M = sum(tot(n)["b0_meet"] for n in NS)
K = sum(tot(n)["b0_cyc"] for n in NS)
MC = sum(tot(n)["bC_meet"] for n in NS)
RS = sum(tot(n)["bC_rescue"] for n in NS)
print("POOLED tapes", T, "meet0", M, "%.4f" % (M / T), "cap0", K, "%.4f" % (K / T),
      "meetC", MC, "%.4f" % (MC / T), "rescued", RS, "%.4f" % (RS / K))

tr = [(r, n, c["bC_trim_max"], c["bC_trim_max"] / c["n"])
      for r, by in R.items() for n, c in by.items() if c.get("bC_trim_max") is not None]
tr.sort(key=lambda x: -x[2])
print("max trims:", tr[0], "| max trims/n:", max(tr, key=lambda x: x[3]))
print("any cap hit under the bit:", any((c.get("bC_cap") or 0) > 0 for by in R.values() for c in by.values()))
print("split violations:", sum((c.get("bC_split_viol") or 0) for by in R.values() for c in by.values()))
