1. What's in the package?

   src/:	the source code of Pestrie and bitmap persistence schemes
   
   test/:   	sample data for testing the code


2. How to use the code?
   
   1. Go to "src" folder and "make";
   2. You can optionally "make install". By default the binaries will be copied to $HOME/bin;
   3. Go to "src/antlr";
   4. Uncompress the file "antlr.ptm.bz2";
   5. Use "pesI antlr.ptm antlr.ptp" to generate the Pestrie persistence file "antlr.ptp". Correspondingly, use "bitI antlr.ptm antlr.ptb" to generate the bitmap persistence file "antlr.ptb";
   6. Use "qtester -t1 antlr.ptp" to test the querying efficiency for alias query. Use "qtester -t1 antlr.ptp basePtrs.in" to calculate the alias pairs with the pointers given in "basePtrs.in". Replacing "antlr.ptp" with "antlr.ptb" will enter the bitmap based querying system;
   7. Type "qtester" without parameters to gain help with other queries.


3. Other usages of this code:

   In fact, Pestrie is an approach to compute, encode, and query the boolean matrix multiplication of A*transpose(A). Pestrie can also be generalized to compute the product A*B for two different matrices. Both are implemented in our code ("pes-self.cc" and "pes-dual.cc").
