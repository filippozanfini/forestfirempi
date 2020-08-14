# Forest Fire with MPI
Forest Fire cellular automata written in C++ with MPI and Allegro5

Forest Fire is defined by four rules which are executed simultaneously:

1. A burning cell turns into an empty cell
2. A tree will burn if at least one neighbor is burning
3. A tree ignites with probability *f* even if no neighbor is burning
4. An empty space fills with a tree with probability *p*.


![Forest Fire](forestfire.gif)
