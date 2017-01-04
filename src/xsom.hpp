#pragma once

#include <xsomConvolution.hpp>
#include <xsomEntity.hpp>
#include <xsomLayer.hpp>
#include <xsomMap.hpp>
#include <xsomNetwork.hpp>
#include <xsomPoint.hpp>
#include <xsomSeq.hpp>
#include <xsomTable.hpp>
#include <xsomUtility.hpp>


/**
 * @example example-000-001-thalamic.cpp  
 * @example example-000-002-cortical.cpp
 * @example example-000-003-both.cpp
 *
 * @example example-001-001-sequencer.cpp  
 *
 * @example example-002-001-tabular.cpp
 * @example example-002-002-convolution.cpp
 * @example example-002-003-gridplot.cpp
 *
 * @example example-003-001-mackeyglass.cpp
 */ 



/**
 * @mainpage
 *
 * @section Overview
 *
 * The xsom library enables the setting up of multi-map kohonen-like
 * architectures. It is based on a c++ functional approach, which may
 * give the code an unusual flavor.
 *
 * @section Fields
 *
 * The main components of xsom are maps. A map can be viewed,
 * mathematically, as a field. It is a function that associates a
 * value to each point of some tolopogical space. For example, an
 * array can be such a function :
 * @code
 * typedef std::array<double,101> Field;
 * ...
 * Field f;
 * a = f[13];
 * c = f[92];
 * d = f[100];
 * e = f[0];
 * @endcode
 *
 * Here, the topological space is [0..99] and the values are
 * scalars. We can also use such an array for coding the contiguous
 * field functions that maps [0,1] to scalar values. The array can
 * code values for f(0/100), f(1/100), f(2/100), ... f(100/100). If
 * one asks for f(.3879), we can interpolate between f(38/100) and
 * f(39/100). In other words, the array f becomes a tabular
 * representation of a function [0,1] -> R, i.e. a function
 * prametrized by all the 101 values in the array.
 *
 * The xsom library provides tools for implementing such tabular
 * functions, but any kind of functions can be used as a field.
 *
 * @section Learning
 *
 * Indeed, in xsom, fields a more than functions. They are indeed
 * adaptive functions. Indeed, if we consider once again tabular-based
 * functions, changing the values in the table enables to change the
 * field function. In a more functional approach, as tabular-based
 * functions are a specific case of parametrized functions, a function
 * is adaptive when the parameters used in its definition can be
 * modified. Such modification can be qualified as "learning".
 *
 * In xsom, learning occurs as follows. One has to define a function
 * "set" (or learn) that takes as argument another function d (our
 * approach is functional). The function "set" should ensure that once
 * it has been called, the field f is modified such as for all
 * positions x, f(x) is close to d(x). This is how this can be
 * implemented for the tabular mapping.
 * 
 * @code
 * typedef unsigned int           Position; // This represents [0..100].
 * typedef double                 Value;
 * typedef std::array<double,101> Field;
 *
 * double desired_at(Position p) {
 *  return sin(.001*p);
 * }
 *
 * template<typename FUN>
 * void set(Field& f, const FUN& d) {
 *   for(Position i = 0; i < 101; ++i)
 *     f[i] = d(i);
 * }
 *
 * Field f;
 * set(f,desired_at); // this is learning.
 * @endcode
 *
 * @section Self-Organizing Maps
 *
 * A self-organizing map (SOM) is indeed a field of prototypes. The
 * topological space used to define the field is the positions in the
 * map (e.g. consider 1D or 2D grids, but also a continuum of such
 * positions). Prototypes can me matched against som input. In xsom,
 * this is how distances between a prototype and the current input is
 * evaluated. First we match each prototype to the current value. Each
 * match produced a scalar, which is hight if the prototype and the
 * input have a good match (i.e. their distance is low). The best
 * matching prototype is the one having the highest match. The result
 * of all matches is a field as well (a field of scalar match values
 * indeed). What is often called the best matching unit (BMU) in
 * SOM-based approaches is here the position where that matching field
 * is maximal.
 *
 * In xsom, multiple SOMs are considered. The consequence is that each
 * computing stage of a SOM process has to be decomposed and that each
 * of such stage can be duplicated in some specific implementation. In
 * xsom, the computing process thus consists in "update" stages and
 * "learn" stages. Let us detail it for the basic SOM computation.
 *
 * @code
 * while(true) {
 *   Sample   xi  = get_new_sample();
 *   Position bmu;
 *   
 *   // ### Update of the prototypes "layer" ###
 *
 *   auto desired_matching_activity_fct = [](Position p) -> Value {
 *     return match(p, prototype(p), xi); // At p, the prototype is matched against the current input.
 *   };
 *   layer_activity_learn(desired_matching_activity_fct);
 *
 *   // ### Update whole Kohonen map ###
 *
 *   map_activity_learn(merge); 
 *   // merge(p) computes a scalar values from all the layer activities at p. 
 *   // This merging is the activity of the map. Here, since we have only one 
 *   // layer, merge(p) = layer_activity(p) is ok.
 *
 *   modify_bmu(bmu, argmax_map_activity()); 
 *   // This changes the bmu (passed by reference) according to the
 *   // position where the map activity (i.e. the merge of all layer
 *   // activities) is the highest. This location is returned by
 *   // argmax_map_activity. As we implement a simple Kohonen map,
 *   // modify_bmu can simply consists in the affectation:
 *   // bmu = argmax_map_activity();
 *   
 * }
 * @endcode
 **/
