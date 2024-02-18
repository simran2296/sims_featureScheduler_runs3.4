

if __name__ == "__main__":
    for day in [2, 3]:
        for filts in ["gz", "griz", "ugrizy"]:
            for rate in [32, 64, 128]:
                print("python too_elab.py --too_days %i --too_rate %i --too_filters %s" % (day, rate, filts))
