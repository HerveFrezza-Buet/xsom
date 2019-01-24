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
 * @example example-002-004-recsom1D.cpp
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
 * value to each element of some tolopogical space. These elements
 * will be referred to as positions in the following. For example, an
 * array can be such a function :
 *
 * @code
 * typedef unsigned int           Position;
 * typedef double                 Value
 * typedef std::array<Value, 101> Field;
 * ...
 * Field f;
 * Position p = 13;
 * Value a = f[p];
 * @endcode
 *
 * Here, the topological space is [0..100] and the values are
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
 * Indeed, in xsom, fields a more than functions. They are actually
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
 * typedef std::array<Value, 101> Field;
 *
 * Value desired_at(Position p) {
 *  return sin(.001*p);
 * }
 *
 * void set_field(Field& f, std::function<Value (Position)> d) {
 *   for(Position i = 0; i < 101; ++i)
 *     f[i] = d(i);
 * }
 *
 * Field f;
 * auto set = std::bind(set_field, std::ref(f), _1);
 * ...
 * set(desired_at); // this is learning... f is handled internally by the set function.
 * @endcode
 * 
 * Using functions as set, which handles internally some field, gives
 * here the flavor of the functions that the user has to provide to
 * xsom.
 *
 * @section Self-Organizing Maps
 *
 * A self-organizing map (SOM) is indeed a field of prototypes. The
 * topological space used to define the field is the positions in the
 * map (e.g. consider 1D or 2D grids, but also a continuum of such
 * positions). Prototypes can me matched against the current input. In xsom,
 * this is how distances between a prototype and the current input is
 * evaluated. First we match each prototype to the current value. Each
 * match produced a scalar, which is high if the prototype and the
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
 * Most functions used in the following embed internal field
 * manipulation, as introduced in previous section.
 *
 * @code
 * Position bmu; // This is what a SOM actually provides as an output.
 *
 * while(true) {
 *   Sample   xi  = get_new_sample();
 *   
 *   // ### Update of the prototypes "layer" ###
 *
 *   auto desired_matching_activity_fct = [](const Position& p) -> Value {
 *     return match(p, prototype(p), xi); 
 *    // At p, the prototype is matched against the current input. A simpler 
 *    // "match(prototype(p), xi);" could have been used, but, for the sake of
 *    // generality, a matching rule that may depend on the location p is 
 *    // rather considered in xsom.
 *   };
 *   layer_activity_learn(desired_matching_activity_fct);
 *
 *   // ### Update whole Kohonen map ###
 *
 *   map_activity_learn(merge); 
 *   // merge(p) computes a scalar values from all the layer activities at p. 
 *   // This merging is the activity of the map that owns the layers. Here, 
 *   // since we have only one layer, using directly layer_activity for the 
 *   // merge function is ok (merge(p) = layer_activity(p)).
 *
 *   modify_bmu(bmu, argmax_map_activity()); 
 *   // This changes the bmu (passed by reference) according to the
 *   // position where the map activity (i.e. the merge of all layer
 *   // activities) is the highest. This location is returned by
 *   // argmax_map_activity. As we implement a simple Kohonen map,
 *   // modify_bmu can simply consists in the affectation:
 *   // bmu = argmax_map_activity();
 *
 *   // ### Learn the prototypes ###
 *
 *   auto new_weight_fct = [](const Position& p) -> Weight {
 *     w = layer_weight(p);
 *     return w + alpha * h(p) * (xi - w);
 *   };
 *   layer_weight_learn(new_weight_fct);
 * }
 * @endcode
 *
 * Indeed, if several weight layers are used, all layers have to be
 * updated in order to compute each layer activity field, then merge combines
 * these activies so that the bmu is updated accordingly. Last,
 * learning is applied to all the layers.
 *
 * The xsom library provides layers and maps. Layers handle weights
 * and matching activities, and the map gathers them and handles the
 * BMU. To sum up, let us recall which functions and data are needed for maps and layers:
 * - Layer elements for computation
 *   - learning rate           (parameter)
 *   - xi                      (data)
 *   - match                   (function)
 *   - h                       (function)
 *   - layer_weight            (function)
 *   - layer_activity_learn    (function)
 *   - layer_weight_learn      (function)
 * - Map  elements for computation
 *   - bmu                  (data)
 *   - map_activity_learn   (function)
 *   - merge                (function)
 *   - argmax_map_activity  (function)
 *   - modify_bmu           (function)
 *
 * The elements listed above are the parameter given to layer and map creation in xsom. 
 **/
