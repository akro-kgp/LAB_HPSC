#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib> // For exit()
#include <cstring> // For strcpy, strcat
#include <ctime>   // For timestamp

using namespace std;

// --- Function Prototypes ---
void gmsh_data_read ( const string& gmsh_filename, int node_dim, int node_num,
  vector<double>& node_x, int element_order, int element_num, vector<int>& element_node );
void gmsh_size_read ( const string& gmsh_filename, int& node_num, int& node_dim,
  int& element_num, int& element_order );
void i4mat_write ( const string& output_filename, int m, int n, const vector<int>& table );
void r8mat_write ( const string& output_filename, int m, int n, const vector<double>& table );
bool s_begin ( const string& s1, const string& s2 );
void timestamp ( );

// --- Main Program ---

int main ( int argc, char *argv[] )
{
  int element_num;
  int element_order;
  string fem_element_filename;
  string fem_node_filename;
  string gmsh_filename;
  int m;
  int node_num;
  string prefix;

  timestamp ( );
  cout << "\n";

  if ( argc <= 1 )
  {
    cout << "\n";
    cout << "  Please enter the filename prefix.\n";
    cin >> prefix;
  }
  else
  {
    prefix = argv[1];
  }

  // --- Create the filenames ---
  gmsh_filename = prefix + ".msh";
  fem_node_filename = prefix + "_nodes.txt";
  fem_element_filename = prefix + "_elements.txt";

  // --- Read GMSH sizes ---
  gmsh_size_read ( gmsh_filename, node_num, m, element_num, element_order );

  // --- Report sizes ---
  cout << "\n";
  cout << "  Size information from GMSH:\n";
  cout << "  Spatial dimension M = " << m << "\n";
  cout << "  Number of nodes NODE_NUM = " << node_num << "\n";
  cout << "  Number of elements ELEMENT_NUM = " << element_num << "\n";
  cout << "  Element order ELEMENT_ORDER = " << element_order << "\n";

  // --- Write info files for the next program ---
  ofstream f1("nodeinfo.txt");
  f1 << node_num << "\n";
  f1 << m << "\n";
  f1.close();

  ofstream f2("eleminfo.txt");
  f2 << element_num << "\n";
  f2 << element_order << "\n";
  f2.close();

  // --- Allocate memory using vectors ---
  vector<double> node_x(m * node_num);
  vector<int> element_node(element_order * element_num);

  // --- Read GMSH data ---
  gmsh_data_read ( gmsh_filename, m, node_num, node_x, element_order,
    element_num, element_node );

  // --- Write FEM data ---
  r8mat_write ( fem_node_filename, m, node_num, node_x );
  i4mat_write ( fem_element_filename, element_order, element_num, element_node );

  // --- Vectors automatically free memory ---

  // --- Terminate ---
  cout << "\n";
  cout << "GMSH_TO_FEM:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}

// --- Function Implementations ---

void gmsh_data_read ( const string& gmsh_filename, int node_dim, int node_num,
  vector<double>& node_x, int element_order, int element_num, vector<int>& element_node )
{
  int dummy;
  int i;
  ifstream input;
  int j;
  int k;
  string line;

  input.open ( gmsh_filename.c_str() );

  if ( !input.is_open() )
  {
    cerr << "\n";
    cerr << "GMSH_DATA_READ - Fatal error!\n";
    cerr << "  Could not open input file \"" << gmsh_filename << "\"\n";
    exit ( 1 );
  }

  // --- Find and read nodes ---
  while ( getline ( input, line ) )
  {
    if ( s_begin ( line, "$Nodes" ) )
    {
      getline ( input, line ); // Read node count line
      j = 0;
      while ( getline ( input, line ) )
      {
        if ( s_begin ( line, "$EndNodes" ) )
        {
          break;
        }
        if ( j < node_num )
        {
            // This is a simplified parser. A more robust one would use stringstream.
            sscanf(line.c_str(), "%d %lf %lf %lf", &dummy, &node_x[0+j*node_dim], &node_x[1+j*node_dim], &node_x[2+j*node_dim]);
            j = j + 1;
        }
      }
      break;
    }
  }

  // --- Find and read elements ---
  while ( getline ( input, line ) )
  {
    if ( s_begin ( line, "$Elements" ) )
    {
      getline ( input, line ); // Read element count line
      j = 0;
      while ( getline ( input, line ) )
      {
        if ( s_begin ( line, "$EndElements" ) )
        {
          break;
        }
        if ( j < element_num )
        {
            // C-style parsing for simplicity and compatibility
            char* line_ptr = (char*)line.c_str();
            int len;
            // Skip the first 5 integers (element number, type, tags)
            for (k=0; k<5; ++k) {
                sscanf(line_ptr, "%d%n", &dummy, &len);
                line_ptr += len;
            }
            // Read the node indices
            for ( i = 0; i < element_order; i++ )
            {
                sscanf(line_ptr, "%d%n", &element_node[i+j*element_order], &len);
                line_ptr += len;
            }
            j = j + 1;
        }
      }
      break;
    }
  }

  input.close ( );

  return;
}

void gmsh_size_read ( const string& gmsh_filename, int& node_num, int& node_dim,
  int& element_num, int& element_order )
{
  int dummy;
  ifstream input;
  string line;
  const double r8_big = 1.0E+30;
  double x, y, z;
  double x_max, x_min, y_max, y_min, z_max, z_min;

  node_num = 0;
  node_dim = 0;
  element_num = 0;
  element_order = 0;

  x_max = - r8_big; x_min = + r8_big;
  y_max = - r8_big; y_min = + r8_big;
  z_max = - r8_big; z_min = + r8_big;

  input.open ( gmsh_filename.c_str() );

  if ( !input.is_open() )
  {
    cerr << "\n";
    cerr << "GMSH_SIZE_READ - Fatal error!\n";
    cerr << "  Could not open input file \"" << gmsh_filename << "\"\n";
    exit ( 1 );
  }

  while ( getline ( input, line ) )
  {
    if ( s_begin ( line, "$Nodes" ) )
    {
      input >> node_num;
      getline ( input, line ); // Consume rest of line

      for(int i = 0; i < node_num; ++i) {
          input >> dummy >> x >> y >> z;
          getline ( input, line ); // Consume rest of line
          if (x < x_min) x_min = x;
          if (x > x_max) x_max = x;
          if (y < y_min) y_min = y;
          if (y > y_max) y_max = y;
          if (z < z_min) z_min = z;
          if (z > z_max) z_max = z;
      }
    }
    else if ( s_begin ( line, "$Elements" ) )
    {
      input >> element_num;
      getline ( input, line ); // Consume rest of line

      if (element_num > 0) {
          getline(input, line); // Read first element line to determine order
          int node_count = 0;
          char* line_ptr = (char*)line.c_str();
          int len;
          // Count integers on the line
          while(sscanf(line_ptr, "%d%n", &dummy, &len) == 1){
              node_count++;
              line_ptr += len;
          }
          element_order = node_count - 5;
      }
      break; // Found what we need
    }
  }

  node_dim = 3;
  if ( z_max == z_min )
  {
    node_dim = 2;
    if ( y_max == y_min )
    {
      node_dim = 1;
    }
  }
  input.close();
}


void i4mat_write ( const string& output_filename, int m, int n, const vector<int>& table )
{
  int i;
  int j;
  ofstream output;

  output.open ( output_filename.c_str() );

  if ( !output.is_open() )
  {
    cerr << "\n";
    cerr << "I4MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the file \"" << output_filename << "\".\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << table[i+j*m];
    }
    output << "\n";
  }

  output.close ( );
}

void r8mat_write ( const string& output_filename, int m, int n, const vector<double>& table )
{
  int i;
  int j;
  ofstream output;

  output.open ( output_filename.c_str() );

  if ( !output.is_open() )
  {
    cerr << "\n";
    cerr << "R8MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file \"" << output_filename << "\".\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << table[i+j*m];
    }
    output << "\n";
  }

  output.close ( );
}

bool s_begin ( const string& s1, const string& s2 )
{
  if (s1.length() < s2.length()) {
      return false;
  }
  return (s1.substr(0, s2.length()) == s2);
}

void timestamp ( )
{
  time_t now = time ( NULL );
  char time_buffer[80];
  strftime(time_buffer, sizeof(time_buffer), "%d %B %Y %I:%M:%S %p", localtime(&now));
  cout << time_buffer << "\n";
}

