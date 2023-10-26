import os
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from math import log
from pandas import read_csv as rcsv
from sys import argv

from utils import user_config


'''Save output as summary results.txt and csv database. Run chart.Charting().'''


# calculate critical temperature of section according to EC3
def temp_crit(coef):
    return 39.19 * log(1 / 0.9674 / coef ** 3.833 - 1) + 482


# calculate rooted mean square error p - probability, n - number of iterations
def rmse(p, n):
    return (p * (1 - p) / n) ** 0.5


def summary(data, t_crit, rset, savepath='.'):
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

    save_list = ['McSteel v0.1.0\n\nResults from {} iterations\n'.format(n_iter)]
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

    with open(os.path.join(savepath, 'results.txt'), 'w') as file:
        file.writelines(save_list)
    print('[OK] Results summary written to TXT file')

    # draw charts
    Plots(data, t_crit, rset, (p_collapse, p_evacfailed)).draw(savepath)

    return err


class Plots:
    def __init__(self, data_frame, t_crit, rset, probs):
        self.results = data_frame
        self.t_crit = t_crit
        self.rset = rset
        self.p_collapse = probs[0]
        self.p_evacfailed = probs[1]
        warnings.filterwarnings('ignore')

    # charts used for risk analysis

    def cdf(self, data, x_crit, y_crit, label, crit_lab):
        sns_plot = sns.distplot(data, hist_kws={'cumulative': True},
                                kde_kws={'cumulative': True, 'label': 'CDF'}, bins=20, axlabel=label)
        plt.axvline(x=x_crit, color='r')
        plt.axhline(y=y_crit, color='r')
        plt.text(x_crit - 0.05 * (plt.axis()[1] - plt.axis()[0]), 0.2, crit_lab, rotation=90)

    def pdf(self, data, x_crit, label, crit_lab):
        sns_plot = sns.distplot(data, kde_kws={'label': 'PDF'}, axlabel=label)
        plt.axvline(x=x_crit, color='r')
        plt.text(x_crit - 0.05 * (plt.axis()[1] - plt.axis()[0]), 0.2, crit_lab, rotation=90)

    def dist(self, savepath='.', type='cdf'):
        try:
            plt.figure(figsize=(12, 4))

            plt.subplot(121)
            if type == 'cdf':
                self.cdf(self.results.temp_max, self.t_crit, 1 - self.p_collapse, 'Temperature [°C]',
                         r'$\theta_{a,cr}$')
            elif type == 'pdf':
                self.pdf(self.results.temp_max, self.t_crit, 'Temperature [°C]', r'$\theta_{a,cr}$')

            plt.subplot(122)
            if type == 'cdf':
                self.cdf(self.results.time_crit[self.results.time_crit > 0], self.rset, self.p_evacfailed,
                         'Time [s]', 'RSET')
            elif type == 'pdf':
                self.pdf(self.results.time_crit[self.results.time_crit > 0], self.rset, 'Time [s]', 'RSET')

        except:
            plt.figure(figsize=(6, 4))
            if type == 'cdf':
                self.cdf(self.results.temp_max, self.t_crit, self.p_collapse, 'Temperature [°C]', r'$\theta_{a,cr}$')
            elif type == 'pdf':
                self.pdf(self.results.temp_max, self.t_crit, 'Temperature [°C]', r'$\theta_{a,cr}$')

        if type == 'cdf':
            savepath = os.path.join(savepath, 'cdf.png')
        elif type == 'pdf':
            savepath = os.path.join(savepath, 'pdf.png')
        plt.savefig(savepath)
        plt.clf()
        plt.close('all')

    def draw(self, path):
        print('Drawing charts...')
        self.dist(path, type='cdf')
        self.dist(path, type='pdf')
        print('[OK] Charts drawn (temp_crit={}, RSET={})'.format(int(self.t_crit), int(self.rset)))

        return 0


'''Drawing densities of fire scenarios and their results'''


class Densities:
    def __init__(self, data_frame, t_crit, rset, probs):
        self.results = data_frame
        self.t_crit = t_crit
        self.rset = rset
        self.p_collapse = probs[0]
        self.p_evacfailed = probs[1]
        warnings.filterwarnings('ignore')

    def basic(self, data, dim=2):
        # basic chart
        pass

    def level_density(self):
        # density in the level of the structure
        pass

    def spatial_density(self):
        # density in given space
        pass

    def console(self):
        # console gui for making density charts
        pass


if __name__ == '__main__':
    user = user_config(argv[1])
    summary(rcsv('{}\\{}_results.csv'.format(user['results_path'], user['case_title'])), temp_crit(user['miu']),
            user['RSET'])
