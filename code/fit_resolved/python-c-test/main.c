#include <python3.8/Python.h>

PyObject* foo(const char* filename)
{
    PyObject* result = PyList_New(0);
    int i;

    for (i = 0; i < 100; ++i)
    {
        PyList_Append(result, PyLong_FromLong(i));
    }

    return result;
}
