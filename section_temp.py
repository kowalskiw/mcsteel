from sys import argv
import numpy as np

# args: path to TEM file, critical temperature
# "python section_temp.py d:\my_sim.gid\hea180.tem 540"


def mean_temp(temfile_path, amb_temp=20):
    nfiber = 0
    is_reading_temp = False
    temp = 0
    t = 0
    section_temp = []

    with open(temfile_path) as file:
        tem = file.readlines()

    for line in tem:
        # set number of elements
        if 'NFIBERBEAM' in line:
            nfiber = int(line.split()[1])

        # set time step value and reset temperature
        elif 'TIME' in line:
            temp = 0
            t = float(line.split()[1])
            is_reading_temp = True

        # try to save previous step mean cross-section temperature
        elif line.startswith('\n'):
            try:
                section_temp.append([t, temp/nfiber])
                is_reading_temp = False
            except UnboundLocalError:
                is_reading_temp = False

        # add temperature in element in certain step
        elif is_reading_temp:
            try:
                fiber_temp = float(line.split()[-1])
                if fiber_temp >= amb_temp:
                    temp += fiber_temp
                else:
                    print('[WARNING] SafirError: Fiber temperature is lower than ambient ({} °C < {} °C)'.format(
                        fiber_temp, amb_temp))
                    raise ChildProcessError
            except (IndexError, UnboundLocalError, ValueError):
                pass

    return np.array(section_temp)


if __name__ == '__main__':
    mtemp_array = mean_temp(argv[1])
    temp_crit = float(argv[2])
    time, temp = zip(*mtemp_array)
    temp_max = max(temp)

    print('_____________________________')
    print('Time, [s] | Temperature, [°C]')
    print('----------|------------------')
    for i in mtemp_array:
        print(int(i[0]), ' '*(10-len(str(i[0]))), '|', round(i[1], 2))

    print('\nMaximum temperature: {}'.format(temp_max))

    if max(temp) > temp_crit:
        for r in mtemp_array:
            if r[1] >= temp_crit:
                crittime = r[0]
                break
        print('\nMean temperature > critical time: {} s'.format(crittime))
    else:
        print('\nCritical temperature has not been reached')
