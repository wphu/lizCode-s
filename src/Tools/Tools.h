/*!
 \page macros Code macros

 The c++ macros used in the code should be placed in the \ref Tools.h file.

 \section caveats Warning, Error and Debug
 All these macros will print to the standard error a tag , the name, line of the source file that caused the call
 as well as the name of the function calling the macro. The arguments contained in the parenthesis will then be
 appended. Arguments can be chained together in c++ stream style (using <tt><<</tt> operator)

 The macro <tt>WARNING("text")</tt> is the most basic and is itended for warnings that should always be present in the code.

 The macro <tt>ERROR("text")</tt> is used to print an error and close the program.

 The macro <tt>DEBUG("text")</tt> can be used in two ways: using just an argument, it will display a debug message
 (similar to <tt>WARNING("text")</tt> ) but it can be used in the form <tt>DEBUG(N,"text")</tt> in this case N is a number and
 represents the debug level starting at which the dubug must be displayed.
 The debug level can be changed int the namelist vie the key <tt>debug</tt>.

 */

#ifndef TOOLS_H
#define TOOLS_H

#include <csignal>
#include <cstdlib>

#include <iostream>


#define __header(__msg,__txt) std::cout << "\t[" << __msg << "] " << __FILE__ << ":" << __LINE__ << " (" \
<< __FUNCTION__ << ") " << __txt << std::endl

#define MESSAGE(__txt)  std::cout << "MESSAGE: ";  std::cout << __txt << std::endl;

#define __PRINTLINE(__num) {MESSAGE(std::string(__num,'-'))}

#define TITLE(...) {MESSAGE(std::endl); MESSAGE(__VA_ARGS__); __PRINTLINE(80);}




#define MESSAGEALL2(__val,__txt) {for (int __i=0;__i<__val;__i++) std::cout << "\t"; MESSAGEALL1(__txt);}
#define MESSAGEALL3(arg1,arg2,arg3,...) arg3
#define MESSAGEALL4(...) MESSAGEALL3(__VA_ARGS__,MESSAGEALL2,MESSAGEALL1,)
#define MESSAGEALL(...) MESSAGEALL4(__VA_ARGS__)(__VA_ARGS__)

#define WARNING(__txt) __header("WARNING", __txt);



#ifdef  __DEBUG

//#warning "DEBUG MODE "
extern int debug_level;

#define DEBUG1(__txt) {if(debug_level>=0) __header("DEBUG", __txt);}
#define DEBUG2(__val,__txt) if(((debug_level<0) && __val==-debug_level) || ((debug_level>=0) && __val<=debug_level)) __header("DEBUG "<<__val, __txt)
#define DEBUG3(arg1,arg2,arg3,...) arg3
#define DEBUG4(...) DEBUG3(__VA_ARGS__,DEBUG2,DEBUG1,)
#define DEBUG(...) DEBUG4(__VA_ARGS__)(__VA_ARGS__)

#define ERROR(__txt) {__header("ERROR", __txt); exit(0);}

#define DEBUGEXEC(...) __VA_ARGS__
#define RELEASEEXEC(...)


#else // __DEBUG

#define DEBUG(...)
#define DEBUGEXEC(...)
#define RELEASEEXEC(...) __VA_ARGS__
#define ERROR(__txt) {__header("ERROR", __txt); exit(0);}

#define HEREIAM(...)

#endif // __DEBUG

class Tools {
 public:
  static void printMemFootPrint(std::string tag);
};





#endif
