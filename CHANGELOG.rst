1.2.0
=============================

 - Avoid forcing the user to know how to use SQLAlchemy
 - Modularize the code, including making all 3D classes subsets of their
   1D classes
 - Add docstrings, ready for Sphinx
 - Allow overlapping associateds, with singlephase() deciding which one to keep
   (for the case where there are a lot of conflicting P-S times)?
      
New methods:
------------

  - associator.assoc1d:
    * count_associated()
    * count_candidates()
    * count_picks()
    * count_stations()
    * get_associated()
    * get_candidates()
    * get_picks()
    * get_stations()
    * print_counts()
    * print_tables()
    * print_table()
    * plot_section().  Like plot1D.section_plot() but don't need to make a
                       new object, also plots all picks (not just associated)

  - associator.tables1d:
    * Pick.from_obspy()
    * Pick.to_obspy()
    * PickModified.to_obspy()

  - associator.misc (new file):
    * isoformat_digits()

  - associator.tt_stations_1d:
    * print_station_table()
    * print_traveltime_table()

  - associator.plot1d: