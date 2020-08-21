#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <sstream>

#include "sampen_calculator.h"

using std::ostringstream;

vector<int> List2Vector(PyObject *py_obj) 
{
    vector<int> result;
    if (PyList_Check(py_obj))
    {
        for (Py_ssize_t i = 0; i < PyList_Size(py_obj); i++)
        {
            PyObject *value = PyList_GetItem(py_obj, i);
            long v = PyLong_AsLong(value);
            if (v == -1) 
            {
                if (PyErr_ExceptionMatches(PyExc_TypeError)) 
                {
                    ostringstream oss;
                    oss << "An element of non integer type has been encountered ";
                    oss << "(index: " << i << "). ";
                    PyErr_SetString(PyExc_TypeError, oss.str().c_str());
                }
            }
            result.push_back(static_cast<int>(v));
        }
    }
    return result;
}

static PyObject * 
sampen_compute_entropy_direct(PyObject *self, PyObject *args)
{
    if (PyTuple_Size(args) != 3) 
    {
        PyErr_SetString(PyExc_TypeError, "Requires exactly three arguments: data, m, r");
        return nullptr;
    }
    vector<int> data = List2Vector(PyTuple_GetItem(args, 0));
    unsigned long m = PyLong_AsUnsignedLong(PyTuple_GetItem(args, 1));
    if (m == static_cast<unsigned long>(-1)) 
    {
        if (PyErr_ExceptionMatches(PyExc_TypeError)) 
        {
            PyErr_SetString(PyExc_TypeError, "m should be a positive integer. ");
            return nullptr;
        }
    }
    if (m > 10) 
    {
        PyErr_SetString(PyExc_TypeError, "m should be an integer ranging from 1 to 10");
        return nullptr;
    }
    long r = PyLong_AsLong(PyTuple_GetItem(args, 2));
    if (r == -1 && PyErr_ExceptionMatches(PyExc_TypeError))
    {
        PyErr_SetString(PyExc_TypeError, "r should be a positive integer. ");
        return nullptr;
    }
    if (r < 0) 
    {
        ostringstream oss;
        oss << "r is supposed to be positive. ";
        PyErr_SetString(PyExc_TypeError, oss.str().c_str());
        return nullptr;
    }
    double a, b;
    double sampen = ComputeSampenDirect(
        data, static_cast<unsigned>(m), static_cast<int>(r), &a, &b);
    
    PyObject *result = PyTuple_New(3);
    if (result == nullptr) 
    {
        PyErr_SetString(PyExc_MemoryError, "Cannot allocate a tuple. ");
        return nullptr;
    }
    PyTuple_SetItem(result, 0, PyFloat_FromDouble(sampen));
    PyTuple_SetItem(result, 1, PyFloat_FromDouble(a));
    PyTuple_SetItem(result, 2, PyFloat_FromDouble(b));
    return result;
}

static PyObject* 
sampen_compute_entropy_qr(PyObject *self, PyObject *args)
{
    if (PyTuple_Size(args) != 6) 
    {
        PyErr_SetString(PyExc_TypeError, "Requires exactly six arguments: data, m, r, presort, sample_size, sample_num");
        return nullptr;
    }
    vector<int> data = List2Vector(PyTuple_GetItem(args, 0));
        unsigned long m = PyLong_AsUnsignedLong(PyTuple_GetItem(args, 1));
    if (m == static_cast<unsigned long>(-1)) 
    {
        if (PyErr_ExceptionMatches(PyExc_TypeError)) 
        {
            PyErr_SetString(PyExc_TypeError, "m should be a positive integer. ");
            return nullptr;
        }
    }
    if (m > 10) 
    {
        PyErr_SetString(PyExc_TypeError, "m should be an integer ranging from 1 to 10");
        return nullptr;
    }
    long r = PyLong_AsLong(PyTuple_GetItem(args, 2));
    if (r == -1 && PyErr_ExceptionMatches(PyExc_TypeError))
    {
        PyErr_SetString(PyExc_TypeError, "r should be a positive integer. ");
        return nullptr;
    }
    if (r < 0) 
    {
        ostringstream oss;
        oss << "r is supposed to be positive. ";
        PyErr_SetString(PyExc_TypeError, oss.str().c_str());
        return nullptr;
    }

    long sample_size = PyLong_AsLong(PyTuple_GetItem(args, 3));
    if (sample_size == -1 && PyErr_ExceptionMatches(PyExc_TypeError))
    {
        PyErr_SetString(PyExc_TypeError, "sample_size should be a positive integer. ");
        return nullptr;
    }
    if (sample_size < 0) 
    {
        ostringstream oss;
        oss << "sample_size (" << sample_size << ") is supposed to be positive. ";
        PyErr_SetString(PyExc_TypeError, oss.str().c_str());
        return nullptr;
    }
    
    long sample_num = PyLong_AsLong(PyTuple_GetItem(args, 4));
    if (sample_num == -1 && PyErr_ExceptionMatches(PyExc_TypeError))
    {
        PyErr_SetString(PyExc_TypeError, "sample_num should be a positive integer. ");
        return nullptr;
    }
    if (sample_num < 0) 
    {
        ostringstream oss;
        oss << "sample_num (" << sample_num << ") is supposed to be positive. ";
        PyErr_SetString(PyExc_TypeError, oss.str().c_str());
        return nullptr;
    }

    int presort = PyObject_IsTrue(PyTuple_GetItem(args, 5));
    if (presort == -1) 
    {
        PyErr_SetString(PyExc_TypeError, "presort should be of bool type. ");
        return nullptr;
    }
    
    double a, b, sampen;
    if (presort) 
    {
        sampen = ComputeSampenQR2(
            data, static_cast<unsigned>(m), static_cast<unsigned>(r), 
            static_cast<unsigned>(sample_size), 
            static_cast<unsigned>(sample_num), &a, &b);
    }
    else 
    {
        sampen = ComputeSampenQR(
            data, static_cast<unsigned>(m), static_cast<unsigned>(r), 
            static_cast<unsigned>(sample_size), 
            static_cast<unsigned>(sample_num), &a, &b);
    }

    PyObject *result = PyTuple_New(3);
    if (result == nullptr) 
    {
        PyErr_SetString(PyExc_MemoryError, "Cannot allocate a tuple. ");
        return nullptr;
    }
    PyTuple_SetItem(result, 0, PyFloat_FromDouble(sampen));
    PyTuple_SetItem(result, 1, PyFloat_FromDouble(a));
    PyTuple_SetItem(result, 2, PyFloat_FromDouble(b));
    return result;    
}

static PyObject* 
sampen_compute_entropy_uniform(PyObject *self, PyObject *args)
{
    if (PyTuple_Size(args) != 5) 
    {
        PyErr_SetString(PyExc_TypeError, "Requires exactly four arguments: data, m, r, presort");
        return nullptr;
    }
    vector<int> data = List2Vector(PyTuple_GetItem(args, 0));
        unsigned long m = PyLong_AsUnsignedLong(PyTuple_GetItem(args, 1));
    if (m == static_cast<unsigned long>(-1)) 
    {
        if (PyErr_ExceptionMatches(PyExc_TypeError)) 
        {
            PyErr_SetString(PyExc_TypeError, "m should be a positive integer. ");
            return nullptr;
        }
    }
    if (m > 10) 
    {
        PyErr_SetString(PyExc_TypeError, "m should be an integer ranging from 1 to 10");
        return nullptr;
    }
    long r = PyLong_AsLong(PyTuple_GetItem(args, 2));
    if (r == -1 && PyErr_ExceptionMatches(PyExc_TypeError))
    {
        PyErr_SetString(PyExc_TypeError, "r should be a positive integer. ");
        return nullptr;
    }
    if (r < 0) 
    {
        ostringstream oss;
        oss << "r is supposed to be positive. ";
        PyErr_SetString(PyExc_TypeError, oss.str().c_str());
        return nullptr;
    }

    long sample_size = PyLong_AsLong(PyTuple_GetItem(args, 3));
    if (sample_size == -1 && PyErr_ExceptionMatches(PyExc_TypeError))
    {
        PyErr_SetString(PyExc_TypeError, "sample_size should be a positive integer. ");
        return nullptr;
    }
    if (sample_size < 0) 
    {
        ostringstream oss;
        oss << "sample_size (" << sample_size << ") is supposed to be positive. ";
        PyErr_SetString(PyExc_TypeError, oss.str().c_str());
        return nullptr;
    }
    
    long sample_num = PyLong_AsLong(PyTuple_GetItem(args, 4));
    if (sample_num == -1 && PyErr_ExceptionMatches(PyExc_TypeError))
    {
        PyErr_SetString(PyExc_TypeError, "sample_num should be a positive integer. ");
        return nullptr;
    }
    if (sample_num < 0) 
    {
        ostringstream oss;
        oss << "sample_num (" << sample_num << ") is supposed to be positive. ";
        PyErr_SetString(PyExc_TypeError, oss.str().c_str());
        return nullptr;
    }

    double a, b;
    double sampen = ComputeSampenUniform(
        data, static_cast<unsigned>(m), static_cast<unsigned>(r), 
        static_cast<unsigned>(sample_size), 
        static_cast<unsigned>(sample_num), &a, &b);

    PyObject *result = PyTuple_New(3);
    if (result == nullptr) 
    {
        PyErr_SetString(PyExc_MemoryError, "Cannot allocate a tuple. ");
        return nullptr;
    }
    PyTuple_SetItem(result, 0, PyFloat_FromDouble(sampen));
    PyTuple_SetItem(result, 1, PyFloat_FromDouble(a));
    PyTuple_SetItem(result, 2, PyFloat_FromDouble(b));
    return result;    
}
static PyMethodDef sampen_methods[] = 
{
    {"compute_sampen_direct", sampen_compute_entropy_direct, METH_VARARGS, "Compute sample entropy using direct method. "}, 
    {"compute_sampen_qr", sampen_compute_entropy_qr, METH_VARARGS, "Compute sample entropy using quasi-Monte Carlo sampling. "}, 
    {"compute_sampen_uniform", sampen_compute_entropy_uniform, METH_VARARGS, "Compute sample entropy using quasi-Monte Carlo sampling. "}, 
    {nullptr, nullptr, 0, nullptr}
};

static PyModuleDef sampen_module =
{
    PyModuleDef_HEAD_INIT, "sampen", "Sample entropy calculation. ", -1, sampen_methods
};

PyMODINIT_FUNC
PyInit_sampen(void) 
{
    return PyModule_Create(&sampen_module);
}
