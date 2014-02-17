/*
  A lightweight test cases management system for experimenting with soot.
  by richardxx, 2012.3
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

using namespace std;

struct benchmark
{
  string name;
  string soot_opts;
  string main_class;
  int n_depend;
  string *pkgs;
  
public:
  void setName( const char* bn )
  {
    name.assign( bn );
  }

  void setSootOpts( const char* opts )
  {
    soot_opts.assign( opts );
  }

  void setMainClass( const char* cl )
  {
    main_class.assign( cl );
  }

  void setNdpkgs( int n )
  {
    n_depend = n;
    pkgs = new string[n];
  }

  void setPackage( const char* pkgName, int i )
  {
    pkgs[i].assign( pkgName );
  }
};

int n_cases;
string global_opts;
string soot_jar;
string out_dir;
string java_opts;
int n_repeat;
string *rep_opts;

const char *arg_list[1024];                   // We permit maximum 1024 parameters
string argBuf[1024];
vector<benchmark*> bm_list;

/*
  The input file format should be:

  n_cases:      the number of test cases
  soot_jar:     the jar file of the soot framework
  out_dir:      the output log files directory
  java_opt:     the options for jvm
  gloal_opt:    the soot options for all test cases

  n_repeat:     the number of repeating times for executing the test cases
  opt_rep1
  opt_rep2
  .....
  opt_repN:     the specific options for each round the whole set test cases are executed

case 1:
  Name:         the name of this test case
  soot_opts:    soot options for this particular test case
  main_class:   the name of the main class
  n_depend:     the number of depending packages:
  pkg1
  pkg2
  .......
  pkgN          the paths to the depending packages

case 2:
  ..........

Notice:
  1. All the soot options can be empty. An empty option list is indicated by a single "." on that line.
  2. A line headed with # is a comment line, we don't process it.
 */
bool parseInput(const char* file_name)
{
  int k;
  string tmp;
  ifstream ifs( file_name, ifstream::in );
  
  if ( ifs.is_open() == false )
    return false;

  // Global settings
  ifs >> n_cases;
  ifs >> soot_jar;
  ifs >> out_dir;
  getline( ifs, tmp );          // We skip the \n character for the current line of input
  getline( ifs, java_opts );
  getline( ifs, global_opts );

  // Repeating settings
  ifs >> n_repeat;
  rep_opts = new string[n_repeat];
  getline( ifs, tmp );
  for ( int i = 0; i < n_repeat; ++i ) {
    getline( ifs, rep_opts[i] );
  }

  // Descriptors for each benchmark
  for ( int i = 0; i < n_cases; ++i ) {
    benchmark *bm = new benchmark();
    
    // Read benchmark name
    ifs >> tmp;
    bm->setName( tmp.c_str() );
    
    // Read specific soot options
    getline( ifs, tmp );
    getline( ifs, tmp );
    bm->setSootOpts( tmp.c_str() );
    
    // Read main class
    ifs >> tmp;
    bm->setMainClass( tmp.c_str() );
    
    // Read the number of depending packages
    ifs >> k;
    bm->setNdpkgs( k );
    
    for ( int j = 0; j < k; ++j ) {
      ifs >> tmp;
      bm->setPackage( tmp.c_str(), j );
    }

    bm_list.push_back( bm );
  }
  
  ifs.close();
  return true;
}

// For debug purpose
void print_args()
{
  int i = 0;

  while ( arg_list[i] != NULL ) {
    cout << arg_list[i] << " ";
    ++i;
  }

  cout << endl;
}

void execute()
{
  int argN;
  int fl_global, fl_bm;
  stringstream ss;

  // We first process the options for java
  argN = 0;
  arg_list[argN++] = "java";
  
  if ( !(java_opts.size() == 1 &&
	 java_opts[0] == '.' ) ) {
    ss.clear();
    ss << java_opts;
    while ( ss >> argBuf[argN] ) {
      arg_list[ argN ] = argBuf[argN].c_str();
      ++argN;
    }
  }

  arg_list[argN++] = "-jar";
  arg_list[argN++] = soot_jar.c_str();


  // Next, we process the global options for soot
  if ( !( global_opts.size() == 1 &&
	  global_opts[0] == '.') ) {
    ss.clear();
    ss << global_opts;
    while ( ss >> argBuf[argN] ) {
      arg_list[ argN ] = argBuf[argN].c_str();
      ++argN;
    }
  }

  // Checkpoint for the global options
  fl_global = argN;

  for ( int i = 0; i < bm_list.size(); ++i ) {
    benchmark *bm = bm_list[i];
    argN = fl_global;                // We override the following options

    // For each benchmark, we first add the output log name
    if ( !( out_dir.size() == 1 &&
	    out_dir[0] == '.' ) ) {
      arg_list[ argN++ ] = "-p";
      arg_list[ argN++ ] = "cg.spark";
      ss.clear();
      ss << "geom-dump-verbose:" << out_dir << "/" << (bm->name);
      ss >> argBuf[argN];
      arg_list[ argN ] = argBuf[argN].c_str();
      ++argN;
    }

    // Then we process the benchmark specific options
    if ( !( bm->soot_opts.size() == 1 &&
	    bm->soot_opts[0] == '.' ) ) {
      ss.clear();
      ss << bm -> soot_opts;
      while ( ss >> argBuf[argN] ) {
	arg_list[ argN ] = argBuf[argN].c_str();
	++argN;
      }
    }

    // Now we process all the depending packages
    if ( bm->n_depend > 0 ) {
      arg_list[ argN++ ] = "-cp";
      
      // We first compute the length of the package list
      int len = bm->n_depend;
      for ( int j = 0; j < bm->n_depend; ++j )
	len += bm->pkgs[j].size();
      argBuf[ argN ].clear();
      argBuf[ argN ].resize( len );
      
      // All the denpending packages are concantenated
      argBuf[argN].assign( bm->pkgs[0] );
      for ( int j = 1; j < bm->n_depend; ++j ) {
	argBuf[ argN ].append( ":" );
	argBuf[ argN ].append( bm->pkgs[j] );
      }
      
      arg_list[ argN ] = argBuf[argN].c_str();
      ++argN;
    }

    // Now we set the main class and the entry class
    arg_list[ argN++ ] = "-main-class";
    arg_list[ argN++ ] = bm->main_class.c_str();
    arg_list[ argN++ ] = bm->main_class.c_str();      // the soot analysis entry class

    // Checkpoint
    fl_bm = argN;

    for ( int k = 0; k < n_repeat; ++k ) {
      argN = fl_bm;

      // We process the repeating options
      if ( !( rep_opts[k].size() == 1 &&
	      rep_opts[k][0] == '.' ) ) {
	ss.clear();
	ss << rep_opts[k];
	while ( ss >> argBuf[argN] ) {
	  arg_list[ argN ] = argBuf[argN].c_str();
	  ++argN;
	}
      }

      arg_list[ argN++ ] = NULL;

      if ( fork() == 0 ) {
	// child
	execv( "/opt/jrockit28.1/bin/java", (char* const*)arg_list );
      }
      else {
	// Father
	int status;
	cout << endl;
	cout << "~~~~~~~~~~~ Processing " << bm->name  << " in round " << k << " ~~~~~~~~~~~~~" << endl;
	print_args();
	cout.flush();
	wait( &status );
      }
    }
  }
}

int main(int argc, char* argv[])
{
  if ( argc < 2 ) {
    cout << "Usage: " << argv[0] << " config_file" << endl;
    return -1;
  }

  if ( parseInput( argv[1] ) == false )
    return -1;
  
  execute();
  return 0;
}
