from cffi import FFI
ffi = FFI()

ffi.cdef('int scale(int, float);')
ffi.set_source('light._distance',
               '''
#include <math.h>

static int scale(int dist, float base)
{
    int result;

    if (base == 1.0){
        return dist;
    }

    if (dist > 0){
        result = (int)(log(dist) / log(base));
        if (result < dist){
           return result;
        }
        else {
           return dist;
        }
    }
    else if (dist < 0){
        result = -1 * (int)(log(-dist) / log(base));
        if (result > dist){
           return result;
        }
        else {
           return dist;
        }
    }
    else {
        return 0;
    }
}
''')

if __name__ == '__main__':
    ffi.compile()
