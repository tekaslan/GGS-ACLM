/**
 * @file
 * @brief Static functions used by geo.c
 * @author Pedro Donato (pdonato@umich.edu) 
*/

/**
 * @brief Functions to Compute Inverse Geodetics per \cite vincenty1975direct .
*/
static void geo_distvincenty(  const struct Pos* const pos1, 
                        const struct Pos* const pos2,                                                            
                        double * const dist, double * const course,
                        const struct GeoOpt * const opt);

/**
 * @brief Functions to Compute Direct Geodetics per \cite vincenty1975direct.
*/
static void geo_nposvincenty(  const struct Pos* const pos1, 
                        struct Pos* const pos2,                                                            
                        const double * const dist, const double * const course,
                        const struct GeoOpt * const opt);

/**
 * @brief Functions to Compute Inverse Geodetics assuming flat earth.
*/
static void geo_distflat(const struct Pos* const pos1, 
                         const struct Pos* const pos2,                                                            
                         double * const dist, double * const course);

/**
 * @brief Functions to Compute Direct Geodetics assuming flat earth.
*/
static void geo_nposflat(  const struct Pos* const pos1, 
                           struct Pos* const pos2,                                                            
                           const double * const dist, 
                           const double * const course);

/**
 * @brief Check is lat and lon values are valid.
*/
static bool geo_checklatlon(const struct Pos* const pos,
                            const struct GeoOpt * const opt);
