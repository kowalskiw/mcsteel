import matplotlib.pyplot as plt
import seaborn as sns
import warnings


class Charting:
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

    def dist(self, type='cdf'):
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
            plt.savefig('dist_p.png')
            plt.clf()
            plt.close('all')
        elif type == 'pdf':
            plt.savefig('dist_d.png')
            plt.clf()
            plt.close('all')

    def draw(self):
        print('Drawing charts...')
        self.dist(type='cdf')
        self.dist(type='pdf')
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
