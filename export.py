from chart import Charting
from math import log
from pandas import read_csv as rcsv
from sys import argv

from fdsafir import user_config


'''Save output as summary results.txt and csv database. Run chart.Charting().'''


# calculate critical temperature of section according to EC3
def temp_crit(coef):
    return 39.19 * log(1 / 0.9674 / coef ** 3.833 - 1) + 482


# calculate rooted mean square error p - probability, n - number of iterations
def rmse(p, n):
    return (p * (1 - p) / n) ** 0.5


def summary(data, t_crit, rset):
    num_nocoll = len(data.time_crit[data.time_crit == 0])
    n_iter = len(data.temp_max)

    def uncertainty(save_list, p, n):
        if save_list[-1][-4:-2] == "0." or save_list[-1][-4:-2] == "1.":
            err = 3 / n
            save_list.append('CI={}\n'.format(err))
        else:
            err = (p * (1 - p) / n) ** 0.5
            save_list.append('RMSE={}\n'.format(err))

        return err, save_list

    print('Saving results...')

    save_list = ['McSteel v0.0.1\n\nResults from {} iterations\n'.format(n_iter)]
    err = [1, 1]  # actual uncertainty of calculation

    # calculating and writing exceeding critical temperature probability and uncertainty to the list
    try:
        p_collapse = len(data.temp_max[data.temp_max >= int(t_crit)]) / len(data.temp_max)
        save_list.append('P(collapse) = {}\n'.format(p_collapse))
    except ZeroDivisionError:
        save_list.append('unable to calculate P(ASET<RSET) and RMSE\n')
        p_collapse = 0
    err[0], save_list = uncertainty(save_list, p_collapse, n_iter)

    # calculating and writing ASET<RSET probability and uncertainty to the list
    try:
        p_evacfailed = (len(data.time_crit[data.time_crit <= int(rset)]) - num_nocoll) / (
                len(data.time_crit) - num_nocoll)
        save_list.append('P(ASET < RSET) = {}\n'.format(p_evacfailed))
    except ZeroDivisionError:
        save_list.append('unable to calculate P(ASET<RSET) and RMSE\n')
        p_evacfailed = 0
    err[1], save_list = uncertainty(save_list, p_evacfailed, n_iter)

    with open('results.txt', 'w') as file:
        file.writelines(save_list)
    print('[OK] Results summary written to TXT file')

    # draw charts
    Charting(data, t_crit, rset, (p_collapse, p_evacfailed)).draw()

    # check if uncertainty is low enough to stop calculations
    if 0 < err[0] < 0.001 and 0 < err[1] < 0.001:
        return True
    else:
        return False


if __name__ == '__main__':
    user = user_config(argv[1])
    summary(rcsv('{}\\{}_results.csv'.format(user['results_path'], user['case_title'])), temp_crit(user['miu']),
            user['RSET'])
