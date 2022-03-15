from sys import argv
from time import time as sec

import mc
import multi
import selsim
from fdsafir2 import out, Config

global outpth

if __name__ == '__main__':
    print(out(outpth, 'McSteel  Copyright (C) 2022  Kowalski W.'
                      '\nThis program comes with ABSOLUTELY NO WARRANTY.'
                      '\nThis is free software, and you are welcome to redistribute it under certain conditions.'
                      '\nSee GPLv3.0 for details (https://www.gnu.org/licenses/gpl-3.0.html).\n'))

    # import configuration files
    cfg = Config(argv[1])

    # run multisimulation of thermal response until stop condition is achieved
    t = sec()
    status = False
    while not status:
        preparations = mc.PrepareMulti(cfg)
        calc_set = preparations.do()

        q = multi.Queue(calc_set)
        status = q.run()

    criteria = ['scenarios number', 'achieving RMSE']
    print(out(outpth, f'[OK] Thermal response multisimulation finished because of {criteria[q.status-1]} limit in'
                      f' {round(sec() - t, 2)} s'))

    # prepare and run the worst scenarios for complex analysis
    t = sec()
    selection = selsim.Selection(cfg)
    calc_set = selection.select()
    q = multi.Queue(selection, mode='full')
    q.run()
    print(out(outpth, f'[OK] Full thermal and mechanical response multisimulation finished in {round(sec() - t, 2)} s'))

    print(out(outpth, 'Thank you for using mcsteel package :)\n'
                      '\nVisit project GitHub site: https://github.com/kowalskiw/mcsteel and contribute!\n'))
