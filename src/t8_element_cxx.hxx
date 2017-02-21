/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** \file t8_element.h
 * This file defines basic operations on an element in a refinement tree.
 *
 * All operations work for all element classes by providing a virtual function table.
 * For each element class, one implementation of the type and virtual table is required.
 */

#ifndef T8_ELEMENT_CXX_HXX
#define T8_ELEMENT_CXX_HXX

#include <sc_refcount.h>
#include <t8_eclass.h>

/** Opaque structure for a generic element, only used as pointer.
 * Implementations are free to cast it to their internal data structure.
 */
typedef struct t8_element t8_element_t;

/** This typedef holds virtual functions for a particular element class. */
struct t8_eclass_scheme_c
{
protected:
  /** This scheme defines the operations for a particular element class. */
  t8_eclass_t eclass;         /**< The element class */
  size_t              element_size; /**< The size in bytes of an element of class \a eclass */
  void               *ts_context;               /**< Anonymous implementation context. */

public:
  /** The destructor. It does nothing but has to be defined since
   * we may want to delete an eclass_scheme that is actually inherited
   * (for example t8_default_scheme_quad) and providing and implementation
   * for the destructor ensures that the
   * destructor of the child class will be executed. */
                      virtual ~ t8_eclass_scheme_c ()
  {
  }

  /** The virtual table for a particular implementation of an element class. */

  /** Return the size of any element of a given class.
   * \return                      The size of an element of class \b ts.
   * We provide a default implementation of this routine that should suffice
   * for most use cases.
   */
  virtual size_t      t8_element_size (void);

  /** Return the maximum allowed level for any element of a given class.
   * \return                      The maximum allowed level for elements of class \b ts.
   */
  virtual int         t8_element_maxlevel (void) = 0;

  /** Return the type of each child in the ordering of the implementation.
   * \param [in] childid  Must be between 0 and the number of children (exclusive).
   *                      The number of children is defined in \a t8_eclass_num_children.
   * \return              The type for the given child.
   */
  virtual t8_eclass_t t8_element_child_eclass (int childid) = 0;

  /** Return the level of a particular element.
   * \param [in] elem    The element whose level should be returned.
   * \return             The level of \b elem.
   */
  virtual int         t8_element_level (const t8_element_t * elem) = 0;

  /** Copy all entries of \b source to \b dest. \b dest must be an existing
   *  element. No memory is allocated by this function.
   * \param [in] source The element whose entries will be copied to \b dest.
   * \param [in,out] dest This element's entries will be overwritted with the
   *                    entries of \b source.
   */
  virtual void        t8_element_copy (const t8_element_t * source,
                                       t8_element_t * dest) = 0;

  /** Compare two elements.
   * \param [in] elem1  The first element.
   * \param [in] elem2  The second element.
   * \return       negativ if elem1 < elem2, zero if elem1 equals elem2
   *               and positiv if elem1 > elem2.
   *  If elem2 is a copy of elem1 then the elements are equal.
   */
  virtual int         t8_element_compare (const t8_element_t * elem1,
                                          const t8_element_t * elem2) = 0;

  /** Compute the parent of a given element \b elem and store it in \b parent.
   *  \b parent needs to be an existing element. No memory is allocated by this function.
   *  \b elem and \b parent can point to the same element, then the entries of
   *  \b elem are overwritten by the ones of its parent.
   * \param [in] elem   The element whose parent will be computed.
   * \param [in,out] parent This element's entries will be overwritten by those
   *                    of \b elem's parent.
   *                    The storage for this element must exist
   *                    and match the element class of the parent.
   *                    For a pyramid, for example, it may be either a
   *                    tetrahedron or a pyramid depending on \b elem's childid.
   */
  virtual void        t8_element_parent (const t8_element_t * elem,
                                         t8_element_t * parent) = 0;

  /** Compute a specific sibling of a given element \b elem and store it in \b sibling.
   *  \b sibling needs to be an existing element. No memory is allocated by this function.
   *  \b elem and \b sibling can point to the same element, then the entries of
   *  \b elem are overwritten by the ones of its i-th sibling.
   * \param [in] elem   The element whose parent will be computed.
   * \param [in] sibid  The id of the sibling computed.
   * \param [in,out] sibling This element's entries will be overwritten by those
   *                    of \b elem's sibid-th sibling.
   *                    The storage for this element must exist
   *                    and match the element class of the sibling.
   *                    For a pyramid, for example, it may be either a
   *                    tetrahedron or a pyramid depending on \b sibid.
   */
  virtual void        t8_element_sibling (const t8_element_t * elem,
                                          int sibid,
                                          t8_element_t * sibling) = 0;

  /** Construct the child element of a given number.
   * \param [in] elem     This must be a valid element, bigger than maxlevel.
   * \param [in] childid  The number of the child to construct.
   * \param [in,out] child        The storage for this element must exist
   *                              and match the element class of the child.
   *                              For a pyramid, for example, it may be either a
   *                              tetrahedron or a pyramid depending on \a childid.
   *                              This can be checked by \a t8_element_child_eclass.
   *                              On output, a valid element.
   * \see t8_element_child_eclass
   */
  virtual void        t8_element_child (const t8_element_t * elem,
                                        int childid, t8_element_t * child) =
    0;

  /** Construct all children of a given element.
   * \param [in] elem     This must be a valid element, bigger than maxlevel.
   * \param [in] length   The length of the output array \a c must match
   *                      the number of children.
   * \param [in,out] c    The storage for these \a length elements must exist
   *                      and match the element class in the children's ordering.
   *                      On output, all children are valid.
   * \see t8_eclass_num_children
   * \see t8_element_child_eclass
   */
  virtual void        t8_element_children (const t8_element_t * elem,
                                           int length, t8_element_t * c[]) =
    0;

  /** Compute the child id of an element.
   * \param [in] elem     This must be a valid element.
   * \return              The child id of elem.
   */
  virtual int         t8_element_child_id (const t8_element_t * elem) = 0;

  /** Query whether a given set of elements is a family or not.
   * \param [in] fam      An array of as many elements as an element of class
   *                      \b ts has children.
   * \return              Zero if \b fam is not a family, nonzero if it is.
   */
  virtual int         t8_element_is_family (t8_element_t ** fam) = 0;

  /* TODO: This could be problematic for pyramids, since elem1 and elem2
   *       could be of different classes. Would need two eclass_schemes as input */
  /** Compute the nearest common ancestor of two elements. That is,
   * the element with highest level that still has both given elements as
   * descendants.
   * \param [in] elem1    The first of the two input elements.
   * \param [in] elem2    The second of the two input elements.
   * \param [in,out] nca  The storage for this element must exist
   *                      and match the element class of the child.
   *                      On output the unique nearest common ancestor of
   *                      \b elem1 and \b elem2.
   */
  virtual void        t8_element_nca (const t8_element_t * elem1,
                                      const t8_element_t * elem2,
                                      t8_element_t * nca) = 0;
  /* TODO: document better */
/** Construct all codimension-one boundary elements of a given element. */
  virtual void        t8_element_boundary (const t8_element_t * elem,
                                           int min_dim, int length,
                                           t8_element_t ** boundary) = 0;

  /** Initialize the entries of an allocated element according to a
   *  given linear id in a uniform refinement.
   * \param [in,out] elem The element whose entries will be set.
   * \param [in] level    The level of the uniform refinement to consider.
   * \param [in] id       The linear id.
   *                      id must fulfil 0 <= id < 'number of leafs in the uniform refinement'
   */
  virtual void        t8_element_set_linear_id (t8_element_t * elem,
                                                int level, uint64_t id) = 0;

  /** Compute the linear id of a given element in a hypothetical uniform
   * refinement of a given level.
   * \param [in] elem     The element whose id we compute.
   * \param [in] level    The level of the uniform refinement to consider.
   * \return              The linear id of the element.
   */
  virtual u_int64_t   t8_element_get_linear_id (const
                                                t8_element_t *
                                                elem, int level) = 0;

  /** Compute the first descendant of a given element.
   * \param [in] elem     The element whose descendant is computed.
   * \param [out] desc    The first element in a uniform refinement of \a elem
   *                      of the maximum possible level.
   */
  virtual void        t8_element_first_descendant (const t8_element_t *
                                                   elem,
                                                   t8_element_t * desc) = 0;

  /** Compute the last descendant of a given element.
   * \param [in] elem     The element whose descendant is computed.
   * \param [out] desc    The last element in a uniform refinement of \a elem
   *                      of the maximum possible level.
   */
  virtual void        t8_element_last_descendant (const t8_element_t *
                                                  elem,
                                                  t8_element_t * desc) = 0;

  /** Construct the successor in a uniform refinement of a given element.
   * \param [in] elem1    The element whose successor should be constructed.
   * \param [in,out] elem2  The element whose entries will be set.
   * \param [in] level    The level of the uniform refinement to consider.
   */
  virtual void        t8_element_successor (const t8_element_t * t,
                                            t8_element_t * s, int level) = 0;

/** Get the integer coordinates of the anchor node of an element */
  /* TODO: better document this */
  virtual void        t8_element_anchor (const t8_element_t * elem,
                                         int anchor[3]) = 0;

  /** Compute the root lenght of a given element, that is the length of
   * its level 0 ancestor.
   * \param [in] elem     The element whose root length should be computed.
   * \return              The root length of \a elem
   */
  virtual int         t8_element_root_len (const t8_element_t * elem) = 0;

  /** Return a pointer to a t8_element in an array indexed by a size_t.
   * \param [in] array    The \ref sc_array storing \t t8_element_t pointers.
   * \param [in] it       The index of the element that should be returned.
   * \return              A pointer to the it-th element in \b array.
   * We provide a default implementation of this routine that should suffice
   * for most use cases.
   */
  virtual t8_element_t *t8_element_array_index (sc_array_t * array,
                                                size_t it);

  /** Allocate memory for an array of elements of a given class.
   * \param [in] length   The number of elements to be allocated.
   * \param [in,out] elems On input an array of \b length many unallocated
   *                      element pointers.
   *                      On output all these pointers will point to an allocated
   *                      and uninitialized element.
   */
  virtual void        t8_element_new (int length, t8_element_t ** elem) = 0;

  /** Deallocate an array of elements.
   * \param [in] ts       The virtual table for this element class.
   * \param [in] length   The number of elements in the array.
   * \param [in,out] elems On input an array of \b length many allocated
   *                      element pointers.
   *                      On output all these pointers will be freed.
   *                      \b elem itself will not be freed by this function.
   */
  virtual void        t8_element_destroy (int length,
                                          t8_element_t ** elem) = 0;
};

/** The scheme holds implementations for one or more element classes. */
struct t8_scheme_cxx
{
  /** Reference counter for this scheme. */
  sc_refcount_t       rc;

  /** This array holds one virtual table per element class. */
  t8_eclass_scheme_c *eclass_schemes[T8_ECLASS_COUNT];
};

/** Increase the reference counter of a scheme.
 * \param [in,out] scheme       On input, this scheme must be alive, that is,
 *                              exist with positive reference count.
 */
void                t8_scheme_cxx_ref (t8_scheme_cxx_t * scheme);

/** Decrease the reference counter of a scheme.
 * If the counter reaches zero, this scheme is destroyed.
 * \param [in,out] pscheme      On input, the scheme pointed to must exist
 *                              with positive reference count.  If the
 *                              reference count reaches zero, the scheme is
 *                              destroyed and this pointer set to NULL.
 *                              Otherwise, the pointer is not changed and
 *                              the scheme is not modified in other ways.
 */
void                t8_scheme_cxx_unref (t8_scheme_cxx_t ** pscheme);

/* TODO: Copy the doxygen comments to the class definition above,
 * then delete all the functions below */
#if 0
/** Destroy an implementation of a particular element class. */
void                t8_eclass_scheme_destroy (t8_eclass_scheme_t * ts);

/* TODO: This function does not exist yet in the eclass_scheme class */
/** Allocate a set of elements suitable for the boundary of a given class.
 * \param [in] scheme           Defines the implementation of the element class.
 * \param [in] theclass         The element class whose boundary we want.
 * \param [in] min_dim          Ignore boundary points of lesser dimension.
 * \param [in] length           Must be equal to the return value
 *                              of \ref t8_eclass_count_boundary.
 * \param [in,out] boundary     On input, array of element pointers of at
 *                              least length \b length.  Filled on output.
 */
void                t8_eclass_boundary_new (t8_scheme_t * scheme,
                                            t8_eclass_t theclass, int min_dim,
                                            int length,
                                            t8_element_t ** boundary);

/* TODO: This function does not exist yet in the eclass_scheme class */
/** Destroy a set of elements suitable for the boundary of a given class.
 * \param [in] scheme           Defines the implementation of the element class.
 * \param [in] theclass         The element class whose boundary we have.
 * \param [in] min_dim          Ignore boundary points of lesser dimension.
 * \param [in] length           Must be equal to the return value
 *                              of \ref t8_eclass_count_boundary.
 * \param [in,out] boundary     Array of element pointers holding elements
 *                              as created by \ref t8_eclass_boundary_new.
 *                              The elements are destroyed by this function.
 */
void                t8_eclass_boundary_destroy (t8_scheme_t * scheme,
                                                t8_eclass_t theclass,
                                                int min_dim, int length,
                                                t8_element_t ** boundary);
#endif /* if 0 */

#endif /* !T8_ELEMENT_CXX_HXX */
