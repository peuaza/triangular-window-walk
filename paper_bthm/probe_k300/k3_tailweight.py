"""The tail's weight: the share of meeting tapes costing more than 100 n and 300 n.

Backs the sentence in bthm_jalc.tex §classes: "the share of meeting tapes costing more than
300 n grows from 0.21 to 0.64", which is why e's local exponent on the band rules (1.5-1.8)
is a crossover and not an order.
"""
import glob
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NS = (512, 1024, 2048, 4096)


def main():
    print(f"{'rule':<10}{'c>100n @512..4096':>26}{'c>300n @512..4096':>26}")
    for rule in sorted({os.path.basename(f).split('_')[1] for f in glob.glob(os.path.join(HERE, 'k3_*.json'))}):
        r100, r300 = [], []
        for n in NS:
            f = os.path.join(HERE, f'k3_{rule}_{n}.json')
            if not os.path.exists(f):
                r100.append(float('nan'))
                r300.append(float('nan'))
                continue
            j = json.load(open(f))
            a = np.array(j['reads'], float)
            r100.append(float(np.mean(a > 100 * n)))
            r300.append(float(np.mean(a > 300 * n)))
        print(f'{rule:<10}{" ".join(f"{x:6.3f}" for x in r100):>26}'
              f'{" ".join(f"{x:6.3f}" for x in r300):>26}')


if __name__ == '__main__':
    main()
