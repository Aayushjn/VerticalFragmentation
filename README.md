# Vertical Fragmentation

The code performs **Vertical Fragmentation** based on a set of inputs provided via the [input.txt](test_1.txt) file.

The file format is as follows:
- attribute count (a)
- query count (q)
- site count (s)
- the next _q_ lines represent the **attribute usage matrix** of size q x a
- the next _q_ lines represent the **query frequency matrix** of size q x s
- the next _q_ lines represent the **query cost matrix** of size q x s