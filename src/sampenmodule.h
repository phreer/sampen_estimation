#ifndef __SAMPENMODULE_H__
#define __SAMPENMODULE_H__

#include <Python.h>

PyObject* sampen_compute_entropy_direct(PyObject *self, PyObject *args);
PyObject* sampen_compute_entropy_qr(PyObject *self, PyObject *args);
PyObject* sampen_compute_entropy_uniform(PyObject *self, PyObject *args);
extern "C"
{
PyMODINIT_FUNC PyInit_sampen(void);
}
#endif // __SAMPENMODULE_H__

