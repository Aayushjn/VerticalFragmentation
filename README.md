# Vertical Fragmentation

The code performs **Vertical Fragmentation** based on a set of inputs provided via the [input.txt](test_1.txt) file.

_Requires numpy and Python 3.8+ to run (uses [Walrus operator](https://realpython.com/lessons/assignment-expressions/))_

_If not using Python 3.8, then modify line 115 in [partition.py](partition.py) to handle variable assignment._

The .txt file format is as follows:
- attribute count (a)
- query count (q)
- site count (s)
- the next _q_ lines represent the **attribute usage matrix** of size q x a
- the next _q_ lines represent the **query frequency matrix** of size q x s
- the next _q_ lines represent the **query cost matrix** of size q x s

Run with `python main.py <file>.txt`