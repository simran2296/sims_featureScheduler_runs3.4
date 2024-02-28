

if __name__ == "__main__":
    for day in [2, 3]:
        for filts in ["gz", "giz", "griz", "ugrizy"]:
            for rate in [128, 256]:
                print("python ../too_elab/too_elab.py --too_days %i --too_rate %i --too_filters %s" % (day, rate, filts))
