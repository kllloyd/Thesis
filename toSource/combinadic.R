combinadic <- function(n, r, i) {
  # FROM:
  #   LSPM: The Leverage Space Portfolio Modeler
  #
  #   Copyright (C) 2009-2010  Soren Macbeth, Joshua Ulrich, and Ralph Vince
  #
  #   This program is free software: you can redistribute it and/or modify
  #   it under the terms of the GNU General Public License as published by
  #   the Free Software Foundation, either version 3 of the License, or
  #   (at your option) any later version.
  #
  #   This program is distributed in the hope that it will be useful,
  #   but WITHOUT ANY WARRANTY; without even the implied warranty of
  #   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  #   GNU General Public License for more details.
  #
  #   You should have received a copy of the GNU General Public License
  #   along with this program.  If not, see <http://www.gnu.org/licenses/>.

  # http://msdn.microsoft.com/en-us/library/aa289166(VS.71).aspx
  # http://en.wikipedia.org/wiki/Combinadic

  if(i < 1 | i > choose(n,r)) stop("'i' must be 0 < i <= n!/(n-r)!")

  largestV <- function(n, r, i) {
    #v <- n-1
    v <- n                                  # Adjusted for one-based indexing
    #while(choose(v,r) > i) v <- v-1
    while(choose(v,r) >= i) v <- v-1        # Adjusted for one-based indexing
    return(v)
  }

  res <- rep(NA,r)
  for(j in 1:r) {
    res[j] <- largestV(n,r,i)
    i <- i-choose(res[j],r)
    n <- res[j]
    r <- r-1
  }
  res <- res + 1

  return(res)
}