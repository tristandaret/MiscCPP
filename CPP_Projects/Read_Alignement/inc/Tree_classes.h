#ifndef MYCLASS_H
#define MYCLASS_H

#include <TROOT.h>
#include <vector>

class ClusterT : public TObject {
public:
   double x;
   double y;
   double z;
   double ex;
   double ey;
   double ez;

   ClusterT() : x(0.0), y(0.0), z(0.0), ex(0.0), ey(0.0), ez(0.0) {}
   virtual ~ClusterT() {}

   ClassDef(ClusterT, 1)
};

class SubtrackT : public TObject {
public:
   int moduleid;
   std::vector<ClusterT> clusters;

   SubtrackT() : moduleid(0) {}
   virtual ~SubtrackT() { clusters.clear(); }

   ClassDef(SubtrackT, 1)
};

class TrackT : public TObject {
public:
   int trackid;
   int n_subtracks;
   std::vector<SubtrackT> subtracks;
   std::vector<double> params;
   std::vector<double> covmat;

   TrackT() : trackid(0), n_subtracks(0) {}
   virtual ~TrackT()
   {
      subtracks.clear();
      params.clear();
      covmat.clear();
   }

   ClassDef(TrackT, 1)
};

class EventT : public TObject {
public:
   int eventid;
   std::vector<TrackT> tracks;

   EventT() : eventid(0) {}
   virtual ~EventT() { tracks.clear(); }

   ClassDef(EventT, 1)
};

class RunT : public TObject {
public:
   std::vector<EventT> events;

   RunT() {}
   virtual ~RunT() { events.clear(); }

   ClassDef(RunT, 1)
};

#endif