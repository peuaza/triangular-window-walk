# -*- coding: utf-8 -*-
"""
汇总分析：对 _dist/<n>.jsonl 逐个 n 给出论文表 2 所需的全部统计量
  N, sd, median, KS/CvM/AD（蒙特卡洛 p 值）, Beta MLE + CI, 端点质量
输出 markdown 表，可直接贴进论文工作稿。
"""
import io, os, json, glob, sys
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DIST = os.path.join(HERE, '_dist')

def load(path):
    xs = []
    for line in io.open(path, encoding='utf-8'):
        line = line.strip()
        if not line:
            continue
        try:
            d = json.loads(line)
        except Exception:
            continue
        if d.get('x') is not None:
            xs.append(d['x'])
    return np.array(xs, dtype=float)

def cvm(u):
    u = np.sort(u); n = len(u); i = np.arange(1, n + 1)
    return 1 / (12 * n) + np.sum((u - (2 * i - 1) / (2 * n)) ** 2)

def ad(u):
    u = np.sort(u); n = len(u); i = np.arange(1, n + 1)
    return -n - np.sum((2 * i - 1) / n * (np.log(u) + np.log(1 - u[::-1])))

def gof(X, B=2000, seed=5):
    u = np.clip(stats.arcsine.cdf(X), 1e-12, 1 - 1e-12)
    t1, t2 = cvm(u), ad(u)
    rng = np.random.default_rng(seed)
    n = len(u)
    n1 = np.empty(B); n2 = np.empty(B)
    for b in range(B):
        s = np.clip(np.sort(rng.random(n)), 1e-12, 1 - 1e-12)
        n1[b] = cvm(s); n2[b] = ad(s)
    return float((n1 >= t1).mean()), float((n2 >= t2).mean())

def analyse(X, B=3000, seed=7):
    N = len(X)
    out = {'N': N, 'sd': float(X.std(ddof=1)), 'median': float(np.median(X))}
    ks = stats.kstest(X, 'arcsine')
    out['ks'] = float(ks.statistic); out['ks_p'] = float(ks.pvalue)
    p1, p2 = gof(X, B=1500)
    out['cvm_p'] = p1; out['ad_p'] = p2
    inter = X[(X > 0) & (X < 1)]
    a, b, _, _ = stats.beta.fit(inter, floc=0, fscale=1)
    out['alpha'] = float(a); out['beta'] = float(b)
    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(200):
        s = rng.choice(inter, len(inter), replace=True)
        try:
            aa, bb, _, _ = stats.beta.fit(s, floc=0, fscale=1)
            bs.append((aa, bb))
        except Exception:
            pass
    if bs:
        bs = np.array(bs)
        out['alpha_ci'] = [float(np.percentile(bs[:, 0], 2.5)), float(np.percentile(bs[:, 0], 97.5))]
        out['beta_ci'] = [float(np.percentile(bs[:, 1], 2.5)), float(np.percentile(bs[:, 1], 97.5))]
    out['mass_lo'] = float((X < 0.05).mean())
    out['mass_hi'] = float((X > 0.95).mean())
    out['q25'] = float(np.quantile(X, 0.25)); out['q75'] = float(np.quantile(X, 0.75))
    return out

def main():
    files = sorted(glob.glob(os.path.join(DIST, '*.jsonl')))
    # 只看主样本 <n>.jsonl，跳过 _rep
    files = [f for f in files if os.path.basename(f).split('.')[0].isdigit()]
    rows = []
    for f in files:
        n = int(os.path.basename(f).split('.')[0])
        X = load(f)
        if len(X) < 100:
            print('skip n=%d (N=%d)' % (n, len(X)))
            continue
        print('分析 n=%d (N=%d) ...' % (n, len(X)), flush=True)
        r = analyse(X); r['n'] = n
        rows.append(r)
    rows.sort(key=lambda r: r['n'])
    print()
    hdr = "| n | N | sd | median | KS | p | CvM p | AD p | α | β | α 95%CI | β 95%CI | P(X<.05) | P(X>.95) |"
    print(hdr)
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        aci = r.get('alpha_ci'); bci = r.get('beta_ci')
        print("| %d | %d | %.4f | %.4f | %.4f | %.3f | %.3f | %.3f | %.3f | %.3f | [%s] | [%s] | %.4f | %.4f |"
              % (r['n'], r['N'], r['sd'], r['median'], r['ks'], r['ks_p'], r['cvm_p'], r['ad_p'],
                 r['alpha'], r['beta'],
                 ('%.3f, %.3f' % tuple(aci)) if aci else '-',
                 ('%.3f, %.3f' % tuple(bci)) if bci else '-',
                 r['mass_lo'], r['mass_hi']))
    print()
    print("理论值：sd 0.3536 ｜ median 0.5 ｜ P(X<.05)=P(X>.95)=0.1436 ｜ α=β=0.5")
    io.open(os.path.join(DIST, 'analysis.json'), 'w', encoding='utf-8').write(
        json.dumps(rows, indent=2, ensure_ascii=False))

if __name__ == '__main__':
    main()
