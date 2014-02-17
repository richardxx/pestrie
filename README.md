1. What's in the package?
   
   src/:           the source code of our algorithms
   pestrie_exp/:   the data used in our experiments
   soot-dump.jar, paddle-dump.jar:  the soot jar packages for producing geometric encoding based and 1H-objsens based points-to information



2. How to reproduce our experimental result in one-step?
   
   Just type "./run" in the command line prompt.
   1. Source code will be compiled and the binaries will be copied to "bin" folder.
   2. Experimental result for each benchmarking subject will be stored in files "build_log.txt" and "query_log.txt" in the subject's own directory. For example, for subject "antlr", its own directory is "pestrie_exp/antlr06". You can find all the directories under "pestrie_exp".



3. How to verify the result in Table 2?
   
   In table 2, we collect the basic characteristics of the subjects. You can find these results in the "build_log.txt" file grep a simple grep. For example, for "antlr", you can run:

   >grep "Input points-to matrix" build_log.txt

   in the command prompt and following line is shown:

   Input points-to matrix: Pointers = 302560, Objects = 76970

   This line gives the number of pointers/objects in original points-to matrix: 302560(76970). Next, run the command:

   >grep "Encoded Points-to matrix"

   The following is shown:

   Encoded Points-to matrix: rows = 88841, columns = 69868, bits = 1413296, mem = 2680kb

   This line gives the number of non-equivalent pointers/objects: 88841(69868). Now you can see, in table 2 of our paper, we have the record for "antlr" is:

   antlr 75.4K 302560 88841 29.4% 76970 69868 90.8%

   Here, 302560 is the number of pointers in input points-to matrix, while 88841 is the non-equivalent number of pointers. You can verify rest of the results in the similar way.

   

4. How to verify the data in Table 7?

   The data relevant to table 7 are given in "query_log.txt". To obtain the "IsAlias" query time for a subject (e.g. "antlr06"), just type in the command:

   >grep "IsAlias querying" query_log.txt

   The results are:

   



5. How to verify the data in Table 8?

   The data shown in table 8 can be found in "build_log.txt" for in every subject's own folder. To check the index size of a subject (e.g. "antlr06"), just type in the command line:

   > grep "index size is" build_log.txt

   There will show two lines similar to following examples:

   The PesTrie index size is : 2551Kb
   The bitmap compressed index size is : 12977Kb

   To checkout the construction time, type:

   > grep "indexing time" build_log.txt

   This also checks out two lines:

   PesTrie indexing time (memory): 540ms (36688Kb) 
   Bitmap indexing time (memory): 749ms (20608Kb)
