/*
 * File: hacker.cc
 *
 * Here are some hack functions to cheat cython compiler
 */

#include <boost/serialization/shared_ptr.hpp>

/*
 * Since cython does not support dereference,  operator*(), *px of shared_ptr
 * (defined by T& operator*() in shared_ptr.hpp) cannot be exposed in cython.
 * Here, we dereference shared_ptr, then assign a new address to px
 */
template <class T>
void assign_deref_shared_ptr(boost::shared_ptr<T>& dest, T *src)
{
    *dest = src;
}

