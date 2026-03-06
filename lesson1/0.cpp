#include <set>
#include <gmsh.h>

int main(int argc, char **argv)
{
  gmsh::initialize();

  gmsh::model::add("circle");

  double radius = 1;
  
  double lc = 1e-2;
  gmsh::model::geo::addPoint(0,0,0, lc, 10);
  gmsh::model::geo::addPoint(0,radius,0, lc, 1);
  gmsh::model::geo::addPoint(radius,0,0, lc, 2);
  gmsh::model::geo::addPoint(0,-radius,0, lc, 3);
  gmsh::model::geo::addPoint(-radius,0,0, lc, 4);
  
  gmsh::model::geo::addCircleArc(1, 10, 2, 1);
  gmsh::model::geo::addCircleArc(2, 10, 3, 2);
  gmsh::model::geo::addCircleArc(3, 10, 4, 3);
  gmsh::model::geo::addCircleArc(4, 10, 1, 4);

  gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);
  gmsh::model::geo::addPlaneSurface({1}, 1);


  gmsh::model::geo::synchronize();

  gmsh::model::mesh::generate(3);

  gmsh::write("circle.msh");

  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();

  gmsh::finalize();

  return 0;
}

