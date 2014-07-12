1. What's in the package?

   src/:	the source code of Pestrie and bitmap persistence schemes
   
   test/:   	sample data for testing the code


2. How to use the code?
   
   1. Go to "src" folder and "make";
   2. You can optionally "make install". By default the binaries will be copied to $HOME/bin;
   3. Go to "src/antlr";
   4. Uncompress the file "antlr.ptm.bz2";
   5. Use "pesI antlr.ptm antlr.ptp" to generate the Pestrie persistence. Replace "pesI" with "bitI" to generate the bitmap persistence;
   6. use "qtester -t1 antlr.ptp" to test the querying efficiency for alias query. Use "qtester -t1 antlr.ptp basePtrs.in" to calculate the alias pairs with the pointers given in "basePtrs.in";
   7. Type "qtester" to get help with other queries.
