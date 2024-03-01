
if __name__ == "__main__":

    # Should be similar to baseline
    print("python ddf_acor.py --ddf_season_frac 0.2 --low_season_fraction 0 --low_season_rate 0")

    ddf_seasons = [0.1, 0.15, 0.2, 0.25, 0.3]

    low_seaons_fractions = [0.15, 0.2, 0.3, 0.4]
    low_season_rates = [0.1, 0.2, 0.5]

    for ddf_season in ddf_seasons:
        for low_seaons_fraction in low_seaons_fractions:
            if low_seaons_fraction > ddf_season:
                for low_season_rate in low_season_rates:
                    print("python ddf_acor.py --ddf_season_frac %.2f --low_season_fraction %.2f --low_season_rate %.2f" % (ddf_season, low_seaons_fraction, low_season_rate))
