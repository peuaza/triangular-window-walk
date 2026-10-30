# -*- coding: utf-8 -*-
"""
run_all.py —— 一键重跑论文里的全部数字

用法:
    python run_all.py                # 全量（慢，约 1–2 小时）
    python run_all.py --quick        # 只跑秒级脚本（验证环境是否装好）
    python run_all.py --only fano.py

每个脚本的输出写到 _results/<script>.out，
耗时写到 _results/timings.csv。
"""
import argparse, os, subprocess, sys, time, io, csv

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, '_results')
os.makedirs(RESULTS, exist_ok=True)

# (脚本, 论文中的位置, 说明, 是否慢)
STEPS = [
    ('wm_exact.py',      'Table 1',        'exhaustive enumeration, n=4..11',        True),
    ('wm12.py',          'Table 1',        'exhaustive enumeration, n=12',           True),
    ('arcsine_conv.py',  'Table 2',        'sd / median / KS vs arcsine',            True),
    ('confirm_c.py',     'Table 3 (A)',    'occurrence rate, convention A',          True),
    ('quarter.py',       'Table 3 (A)',    'occurrence rate, small n',               True),
    ('fallback_big.py',  'Table 3 (B)',    'occurrence rate, convention B',          True),
    ('p_cond.py',        'Sec. 3.4',       'conditional law of X given ac',          True),
    ('ac_corr3.py',      'Conj. 4',        'correlation between consecutive windows', True),
    ('switch_m.py',      'Table 5',        'reversal rate vs m',                     True),
    ('ac_path.py',       'Sec. 4 Link 1',  'path of the window autocorrelation',     False),
    ('cluster_src.py',   'Sec. 4 Link 2',  'reversal rate vs |ac|',                  False),
    ('excursion.py',     'Sec. 4 Link 3',  'excursions of |ac| below threshold',     True),
    ('fano.py',          'Sec. 4 Link 3',  'Fano factors',                           False),
    ('by_switch.py',     'Sec. 4 Link 4',  'run structure by switching',             True),
    ('maxrun_by_k.py',   'Sec. 4 Link 4',  'run structure',                          False),
    ('burst_pos.py',     'Sec. 4 Link 5',  'burst positions',                        False),
    ('segments.py',      'Sec. 2.6',       'decomposition into segments',            True),
    ('qseg.py',          'Sec. 2.6',       'segment statistics',                     True),
]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true', help='跳过慢脚本')
    ap.add_argument('--only', default=None)
    args = ap.parse_args()

    steps = STEPS
    if args.only:
        steps = [s for s in STEPS if s[0] == args.only]
        if not steps:
            sys.exit('unknown script: ' + args.only)
    elif args.quick:
        steps = [s for s in STEPS if not s[3]]

    rows = []
    t_all = time.time()
    for script, where, desc, slow in steps:
        path = os.path.join(HERE, script)
        if not os.path.exists(path):
            print('[MISSING] %s' % script, flush=True)
            continue
        print('== %-18s %-14s %s' % (script, where, desc), flush=True)
        t0 = time.time()
        with io.open(os.path.join(RESULTS, script + '.out'), 'w', encoding='utf-8') as f:
            p = subprocess.run([sys.executable, '-u', script], cwd=HERE,
                               stdout=f, stderr=subprocess.STDOUT)
        dt = time.time() - t0
        status = 'ok' if p.returncode == 0 else 'exit=%d' % p.returncode
        print('   %s  %.1fs' % (status, dt), flush=True)
        rows.append([script, where, desc, '%.1f' % dt, status])

    with io.open(os.path.join(RESULTS, 'timings.csv'), 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f)
        w.writerow(['script', 'paper_location', 'description', 'seconds', 'status'])
        w.writerows(rows)
    print('\n总耗时 %.1f 分钟，结果在 _results/' % ((time.time() - t_all) / 60))

if __name__ == '__main__':
    main()
