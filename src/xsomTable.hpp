#pragma once

#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>

namespace xsom {

  namespace tab {

    template<typename INDEX, typename POSITION>
    class Mapping {
    public:
      
      POSITION min,max;
      INDEX size;

      typedef INDEX    index_type;
      typedef POSITION position_type;

      Mapping(POSITION pos_min, 
	      POSITION pos_max,
	      INDEX    size) 
	: min(pos_min), max(pos_max) {}
      virtual ~Mapping() {}

      INDEX        pos2index(POSITION      pos)  const {}
      POSITION     index2pos(INDEX         idx)  const {}
      unsigned int pos2rank (POSITION      pos)  const {}
      POSITION     rank2pos (unsigned int rank)  const {}
    };

    template<typename CONTENT, typename MAPPING>
    class Table {
    public:
      
      const MAPPING& mapping;
      std::vector<CONTENT> content;

      Table(const MAPPING& m) 
	: mapping(m), 
	  content(m.size)  {}

      CONTENT operator()(typename MAPPING::position_type pos) const {
	return content[mapping.pos2rank(pos)];
      }
      
      void clear() {
	std::fill(content.begin(), content.end(), 0);
      }

      void learn(std::function<CONTENT (typename MAPPING::position_type)> fct_value_at) {
	unsigned int rank;
	unsigned int size  = mapping.size;
	auto c   = content.begin();
	for(rank = 0; rank < size; ++rank, ++c) {
	  auto pos = mapping.rank2pos(rank);
	  *c = fct_value_at(pos);
	}
      }
      
      typename MAPPING::position_type bmu() const {
	auto it = std::max_element(content.begin(), content.end());
	return mapping.rank2pos((unsigned int)(std::distance(content.begin(), std::max_element(content.begin(), content.end())));
      }
	  
      void write(std::ostream& os) const {
	for(auto& val : content) os << val << '\n';
      }

      void read(std::istream& is) {
	char sep;
	for(auto& val : content) {
	  is >> val;
	  is.get(sep);
	}
      }
      
    };
  }


  
  namespace tab1d {
    

    class Mapping : public xsom::tab::Mapping<unsigned int,double> {
    public:
      typedef unsigned int index_type;
      typedef double       position_type;
      double coefx,coefw;

      Mapping(position_type pos_min, position_type pos_max, unsigned int nb_pos) 
	: xsom::tab::Mapping<unsigned int,double>(pos_min,pos_max,nb_pos) {
	coefw = (nb_pos-1.)/(max-min);
	coefx = (max-min)/(nb_pos-1.);
      }
      
      unsigned int pos2index(double pos)  const override {
	return (unsigned int)((pos-min)*coefw+.5);
      }
      
      double index2pos(unsigned int idx)  const override {
	return idx*coefx+min;
      }
      
      unsigned int pos2rank (double pos)  const override {
	return (unsigned int)((pos-min)*coefw+.5);
      }
      
      double rank2pos (unsigned int rank)  const {
	return rank*coefx+min;
      }
    };

    Mapping mapping(double pos_min, double pos_max, unsigned int nb_pos) {
      return Mapping(pos_min,pos_max,nb_pos);
    }


    template<typename CONTENT>
    class Table : public xsom::tab::Table<CONTENT,Mapping> {
    public:
      const Mapping& mapping;
      std::vector<CONTENT> content;

      using xsom::tab::Table<CONTENT,Mapping>::Table;

      void fill_line(std::vector<ccmpl::Point>& points) {
	
	points.clear();

	unsigned int rank;
	unsigned int size  = mapping.size;
	auto c   = l.content.begin();
	for(idx = 0; idx < size; ++idx, ++c) {
	  pos = l.mapping.rank2pos(idx);
	  points.push_back({pos, *c});
	}
      }
    };

    template<typename CONTENT>
    Table<CONTENT> table(const Mapping& m) {
      return Table<CONTENT>(m);
    }
  }

    
  namespace tab2d {
    class Mapping {
    public:

      typedef xsom::Point2D<double>               position_type;
      typedef xsom::Index2D                       index_type;
      position_type min,max;
      index_type size;

      Mapping(const position_type& pos_min, 
	      const position_type& pos_max) 
	: min(pos_min), max(pos_max) {}
      virtual ~Mapping() {}

      virtual index_type    pos2index(const position_type& pos) const  = 0;
      virtual unsigned int  pos2rank(const position_type& pos)  const  = 0;
      virtual position_type index2pos(const index_type& idx)    const  = 0;
    };

    template<int NB_W, int NB_H>
    class Mapping_ : public Mapping {
    public:
      double coefx,coefy,coefw,coefh;
      
      virtual ~Mapping_() {}

      Mapping_(const position_type& pos_min, const position_type& pos_max) 
	: Mapping(pos_min,pos_max) {
	this->size.w = NB_W;
	this->size.h = NB_H;
	coefw = (NB_W-1.)/(this->max.x-this->min.x);
	coefh = (NB_H-1.)/(this->max.y-this->min.y);
	coefx = (this->max.x-this->min.x)/(NB_W-1.);
	coefy = (this->max.y-this->min.y)/(NB_H-1.);
      }

      virtual index_type pos2index(const position_type& pos) const {
	return {
	  (unsigned int)((pos.x-this->min.x)*coefw+.5),
	    (unsigned int)((pos.y-this->min.y)*coefh+.5)
	    };
      }

      virtual unsigned int pos2rank(const position_type& pos) const {
	return (unsigned int)((pos.y-this->min.y)*coefh+.5)*NB_W+(unsigned int)((pos.x-this->min.x)*coefw+.5);
      }

      virtual position_type index2pos(const index_type& idx) const {
	return {
	  idx.w*coefx+this->min.x,
	    idx.h*coefy+this->min.y,
	    };
      }
    };

    template<int NB_W, int NB_H>
    Mapping_<NB_W,NB_H> mapping(const xsom::Point2D<double>& pos_min,
				const xsom::Point2D<double>& pos_max) {
      return Mapping_<NB_W,NB_H>(pos_min,pos_max);
    }


    template<typename CONTENT>
    class Table {
    public:
      const Mapping& mapping;
      std::vector<CONTENT> content;
      std::function<bool (const xsom::Point2D<double>&)> pos_is_valid;
      xsom::Point2D<double> pos_color_min;
      xsom::Point2D<double> pos_color_max;

      template<typename fctPOS_IS_VALID>
      Table(const Mapping& m,
	    const fctPOS_IS_VALID& fct_pos_is_valid) 
	: mapping(m), 
	  content(m.size.w * m.size.h),
	  pos_is_valid(fct_pos_is_valid),
	  pos_color_min(-1,-1),
	  pos_color_max( 1, 1) {}
      virtual ~Table() {}

      CONTENT operator()(const xsom::Point2D<double>& pos) const {
	return content[mapping.pos2rank(pos)];
      }
      
      void clear() {
	std::fill(content.begin(), content.end(), 0);
      }

      virtual void learn(std::function<CONTENT (const xsom::Point2D<double>&)> fct_value_at) {
	xsom::Index2D idx;
	unsigned int width  = mapping.size.w;
	unsigned int height = mapping.size.h;
	auto c   = content.begin();
	for(idx.h = 0; idx.h < height; ++idx.h)
	  for(idx.w = 0; idx.w < width; ++idx.w, ++c) {
	    xsom::Point2D<double> pos = mapping.index2pos(idx);
	    if(pos_is_valid(pos))
	      *c = fct_value_at(pos);
	  }
      }    
     
      // bmu is a template in order to avoid its systematic
      // rewriting. Indeed, is VALUE does not handle operator<, the
      // bmu code fails.
      template<typename VALUE=CONTENT>
      xsom::Point2D<double> bmu() const {
	unsigned int width  = mapping.size.w;
	unsigned int height = mapping.size.h;
	unsigned int w,h;
	auto it = content.begin();
	unsigned int wmax = 0;
	unsigned int hmax = 0;
	VALUE        max  = *it;
	VALUE        m;
	for(h=0;h<height;++h)
	  for(w=0;w<width;++w,++it) 
	    if((m = *it) > max) {
	      wmax = w;
	      hmax = h;
	      max = m;
	    }
	return mapping.index2pos({wmax,hmax});
      }

      void write(std::ostream& os) const {
	for(auto& val : content) os << val << ' ';
      }

      void read(std::istream& is) {
	char sep;
	for(auto& val : content) {
	  is >> val;
	  is.get(sep);
	}
      }
    };

    template<typename CONTENT,
	     typename fctPOS_IS_VALID>
    Table<CONTENT> table(const Mapping& m,
			 const fctPOS_IS_VALID& fct_pos_is_valid) {
      return Table<CONTENT>(m,fct_pos_is_valid);
    }

    void fill_plot_values(const Table<double>& l, std::vector<ccmpl::ValueAt>& points) {
	
      points.clear();

      xsom::Index2D idx;
      unsigned int width  = l.mapping.size.w;
      unsigned int height = l.mapping.size.h;
      auto c   = l.content.begin();
      for(idx.h = 0; idx.h < height; ++idx.h)
	for(idx.w = 0; idx.w < width; ++idx.w, ++c) {
	  xsom::Point2D<double> pos = l.mapping.index2pos(idx);
	  if(l.pos_is_valid(pos))
	    points.push_back({pos.x,pos.y,*c});
	}
    }


    template<typename CONTENT,typename PointOf>
    void fill_plot_color(const Table<CONTENT>& l, 
			 const PointOf& point_of,
			 std::vector<ccmpl::ColorAt>& points) {
      points.clear();

      xsom::Index2D idx;
      unsigned int width  = l.mapping.size.w;
      unsigned int height = l.mapping.size.h;
      auto c   = l.content.begin();
      for(idx.h = 0; idx.h < height; ++idx.h)
	for(idx.w = 0; idx.w < width; ++idx.w, ++c) {
	  xsom::Point2D<double> pos = l.mapping.index2pos(idx);
	  if(l.pos_is_valid(pos))
	    points.push_back({pos.x,pos.y,ccmpl::color::from_point<double>(point_of(*c),
									   l.pos_color_min.x, 
									   l.pos_color_max.x, 
									   l.pos_color_min.y, 
									   l.pos_color_max.y)});
	}
    }
  }
    
    
}
