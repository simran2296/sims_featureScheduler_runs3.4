

if __name__ == "__main__":

    times = [27, 30, 38, 45, 60]
    scales = [0.9, 1, 1.1, 1.2, 1.5]

    for time in times:
        for scale in scales:
            print('python internal_u.py --u_expt %i --u_n_scale %.1f' % (time, scale))
