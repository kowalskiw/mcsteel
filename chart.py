import matplotlib.pyplot as plt
import seaborn as sns


class Charting:
    def __init__(self, data_frame, t_crit, rset, probs):
        self.results = data_frame
        self.t_crit = t_crit
        self.rset = rset
        self.p_coll = probs[0]
        self.p_evac = probs[1]

    # charts used for risk analysis

    def cdf(self, data, x_crit, y_crit, label, crit_lab):
        sns_plot = sns.distplot(data, hist_kws={'cumulative': True},
                                kde_kws={'cumulative': True, 'label': 'CDF'}, bins=20, axlabel=label)
        plt.axvline(x=x_crit, color='r')
        plt.axhline(y=y_crit, color='r')
        plt.text(x_crit - 0.05 * (plt.axis()[1]-plt.axis()[0]), 0.2, crit_lab, rotation=90)

    def pdf(self, data, x_crit, label, crit_lab):
        sns_plot = sns.distplot(data, kde_kws={'label': 'PDF'}, axlabel=label)
        plt.axvline(x=x_crit, color='r')
        plt.text(x_crit - 0.05 * (plt.axis()[1] - plt.axis()[0]), 0.2, crit_lab, rotation=90)

    def dist(self, type='cdf'):
        try:
            plt.figure(figsize=(12, 4))

            plt.subplot(121)
            if type == 'cdf':
                self.cdf(self.results.temp_max, self.t_crit, self.p_coll, 'Temperature [°C]', r'$\theta_{a,cr}$')
            elif type == 'pdf':
                self.pdf(self.results.temp_max, self.t_crit, 'Temperature [°C]', r'$\theta_{a,cr}$')

            plt.subplot(122)
            if type == 'cdf':
                self.cdf(self.results.time_crit[self.results.time_crit > 0], self.rset, self.p_evac, 'Time [s]', 'RSET')
            elif type == 'pdf':
                self.pdf(self.results.time_crit[self.results.time_crit > 0], self.rset, 'Time [s]', 'RSET')

        except:
            plt.figure(figsize=(6, 4))
            if type == 'cdf':
                self.cdf(self.results.temp_max, self.t_crit, self.p_coll, 'Temperature [°C]', r'$\theta_{a,cr}$')
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
        print(self.results)
        self.dist(type='cdf')
        self.dist(type='pdf')

        return 0
