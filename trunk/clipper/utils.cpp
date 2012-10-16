#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "utils.h"

std::string string_format(const std::string &fmt, ...) {
       int n, size=100;
       std::string str;
       va_list ap;
       while (1) {
       str.resize(size);
       va_start(ap, fmt);
       int n = vsnprintf((char *)str.c_str(), size, fmt.c_str(), ap);
       va_end(ap);
       if (n > -1 && n < size)
           return str;
       if (n > -1)
           size=n+1;
       else
           size*=2;
       }
}

std::vector<char *> split(char* str,const char* delim)
{
    char* token = strtok(str,delim);
    std::vector<char *> result;
    while(token != NULL)
    {
        result.push_back(token);
        token = strtok(NULL,delim);
    }
    return result;
}

