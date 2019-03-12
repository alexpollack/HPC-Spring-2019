# include <cstdlib>
# include <iostream>

using namespace std;

int main ( );
void f ( int n );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TEST01.
//  g++ -O3 val_test01.cpp , then ./a.out , then valgrind --tool=memcheck ./a.out then g++ -g inner-mem.cpp
//  Discussion:
//
//    TEST01 calls F, which has a memory "leak".  This memory leak can be
//    detected by VALGRID.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2011
//
{
  int n = 10;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  C++ version.\n";
  cout << "  A sample code for analysis by VALGRIND.\n";

  f ( n );
//
//  Terminate.
//
  cout << "\n";
  cout << "TEST01\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
//****************************************************************************80

void f ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    F computes N+1 entries of the Fibonacci sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2011
//
{
  int i;
  int *x;

  x = ( int * ) malloc ( n * sizeof ( int ) ); //address is 0 byte after malloc size 40

  x[0] = 1;
  cout << "  " << 0 << "  " << x[0] << "\n";

  x[1] = 1;
  cout << "  " << 1 << "  " << x[1] << "\n";

  for ( i = 2; i < n; i++ ) // **WAS <= n**, changed to be < n not <= n, went out of bounds -> bad access
  {
    x[i] = x[i-1] + x[i-2]; //invalid write of size 4, by call at line 40 because of bad acess
    cout << "  " << i << "  " << x[i] << "\n"; //invalide read of size 4, line 40 because of bad access
  }

    free (x);//delete [] x; //mismatched free/delete, think it needs free because of malloc not new

  return;
}
