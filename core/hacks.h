/*
 * File: hacks.h
 *
 * Some hack functions to cheat cython compiler and Block global variables
 */

#include <stdio.h>
#include <boost/shared_ptr.hpp>

void initialize_default_dmrginp(char *fcidump, std::string& prefix, std::string& sym);
int get_last_site_id();

/*
 * Since cython does not support dereference,  operator*(), *px of shared_ptr
 * (defined by T& operator*() in shared_ptr.hpp) cannot be exposed in cython.
 * Here, we dereference shared_ptr, then assign a new address to px
 */
template <class T>
void assign_deref_shared_ptr(boost::shared_ptr<T>& dest, T& src)
{
    *dest = src;
}
