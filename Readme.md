# AFRA — Alignment-Free Support Values

This program computes support values from distance matrices without the
need for a multiple sequence-alignment.

## Compilation Instructions

Execute the following steps to create your own version of `afra`.

    % autoreconf -i
    % ./configure
    % make
    % make install

You should now have a `afra` executable ready for usage.

    % ./afra --mode quartet foo.mat
    (((Seq0:-0.113017,Seq4:0.645317)100:0.084738,Seq1:-0.069237)50:0.053837,Seq3:0.354563,Seq2:0.093537);

## License

Copyright (C) 2015 - 2016  Fabian Klötzl

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
