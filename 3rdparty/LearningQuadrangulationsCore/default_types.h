#ifndef DEFAULT_TYPES_H
#define DEFAULT_TYPES_H

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

/**
 * @brief small_int_type Integer type used for storing small signed integers - [-2^15+1, 2^15-1].
 */
typedef                short int    small_int_type;

/**
 * @brief u_small_int_type Integer type used for storing small unsigned integers - [0, 2^16-1].
 */
typedef   unsigned     short int    u_small_int_type;

#if __cplusplus >= 201103L || _MSC_VER >= 1800
/**
 * @brief big_int_type Integer type used to store big signed integers - [-2^63+1, 2^63-1] (c++11).
 */
typedef            long long int    big_int_type;
#else
/**
 * @brief big_int_type Integer type used to store big signed integers - [-2^31+1, 2^31-1] (c++98).
 */
typedef                 long int    big_int_type;
#endif

#if __cplusplus >= 201103L || _MSC_VER >= 1800
/**
 * @brief u_big_int_type Integer type used to store big unsigned integers - [0, 2^64-1] (c++11).
 */
typedef   unsigned long long int    u_big_int_type;
#else
/**
 * @brief u_big_int_type Integer type used to store big unsigned integers - [0, 2^32-1] (c++98).
 */
typedef   unsigned      long int    u_big_int_type;
#endif

/**
 * @brief n_corners_type Type used to store the number of corners.
 */
typedef   u_small_int_type    n_corners_type;

/**
 * @brief var_label_type Type used to store the index of the variable into its label.
 */
typedef   u_small_int_type    var_label_type;

/**
 * @brief n_corners_type Type used to store the valence of a vertex.
 */
typedef   u_small_int_type    valence_type;

/**
 * @brief num_type Type used to store generic numbers.
 */
typedef   big_int_type        num_type;

/**
 * @brief count_type Type used to store generic numbers.
 */
typedef   u_big_int_type      count_type;

/**
 * @brief SEPARATOR Character used to separate sides' information inside strings.
 */
static const char   SEPARATOR     = '-';

/**
 * @brief BORDER Character used to identify a border adjacency.
 */
static const char   BORDER        = '#';

/**
 * @brief VARLABEL Character used to label a variable.
 */
static const char   VARLABEL      = 'x';

/**
 * @brief VARSEPARATOR Character used to separate variable's name from variable's value.
 */
static const char   VARSEPARATOR  = '=';

/**
 * @brief VALENCELABEL Character used to refer to the valence type of a vertex.
 */
static const char   VALENCELABEL  = 'v';

/**
 * @brief VALENCESEP Character used to separate the valence type from the occurrences.
 */
static const char   VALENCESEP    = ':';

/**
 * @brief POSFACELABEL Character used to label a face index of a Pos.
 */
static const char   POSFACELABEL  = 'F';

/**
 * @brief POSEDGELABEL Character used to label an edge index of a Pos.
 */
static const char   POSEDGELABEL  = 'E';

/**
 * @brief POSVERTLABEL Character used to label a vertex index of a Pos.
 */
static const char   POSVERTLABEL  = 'V';

/**
 * @brief LABLINDEXSEP Character used to separate an index in a Pos from its label.
 */
static const char   LABLINDEXSEP  = '=';

/**
 * @brief INDSEPARATOR Character used to separate several indices of a Pos.
 */
static const char   INDSEPARATOR  = ';';

/**
 * @brief POSSEPARATOR Character used to separate different Poses.
 */
static const char   POSSEPARATOR  = '-';

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // DEFAULT_TYPES_H
