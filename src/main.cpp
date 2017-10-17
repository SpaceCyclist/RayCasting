#include <Geometry/Texture.h>
#include <Math/Vectorf.h>
#include <Geometry/Ray.h>
#include <Geometry/Triangle.h>
#include <Geometry/CastedRay.h>
#include <stdlib.h>
#include <iostream>
#include <Geometry/RGBColor.h>
#include <Geometry/Material.h>
#include <Geometry/PointLight.h>
#include <Geometry/Camera.h>
#include <Geometry/Cube.h>
#include <Geometry/Disk.h>
#include <Geometry/Cylinder.h>
#include <Geometry/Cone.h>
#include <Visualizer/Visualizer.h>
#include <Geometry/Scene.h>
#include <Geometry/Cornel.h>
#include <Geometry/Loader3ds.h>
#include <Geometry/BoundingBox.h>
#include <omp.h>

/// <summary>
/// The directory of the 3D objetcs
/// </summary>
const std::string m_modelDirectory = "../Models";

using Geometry::RGBColor ;

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	void createGround(Geometry::Scene & scene)
///
/// \brief	Adds a ground to the scene.
///
/// \author	F. Lamarche, Universit� de Rennes 1
/// \date	03/12/2013
///
/// \param [in,out]	scene	The scene.
////////////////////////////////////////////////////////////////////////////////////////////////////
void createGround(Geometry::Scene & scene)
{
	Geometry::BoundingBox sb = scene.getBoundingBox();
	// Non emissive
	Geometry::Material * material = new Geometry::Material(RGBColor(), RGBColor(0.5, 0.5, 0.5), RGBColor(0.5, 0.5, 0.5)*8, 1000.0f); // Non existing material...
	//Geometry::Material * material = new Geometry::Material(RGBColor(), RGBColor(0., 0., 0.), RGBColor(1., 1., 1.), 10000.0f);
	//Geometry::Material * material = new Geometry::Material(RGBColor(), RGBColor(1.0f, 1.0f, 1.0f), RGBColor(0.f, 0.f, 0.f), 100.0f);

	//Geometry::Material * material = new Geometry::Material(RGBColor(), RGBColor(1.0,1.0,1.0), RGBColor(), 1000.0f, RGBColor(0.5, 0.5, 0.5)*200); // Non existing material...

	Geometry::Square square(material);
	Math::Vector3f scaleV = (sb.max() - sb.min()) ;
	double scale = ::std::max(scaleV[0], scaleV[1])*2.0;
	square.scaleX(scale);
	square.scaleY(scale);
	Math::Vector3f center = (sb.min() + sb.max()) / 2.0;
	center[2] = sb.min()[2];
	square.translate(center);
	scene.add(square);
	::std::cout << "Bounding box: " << sb.min() << "/ " << sb.max() << ", scale: "<<scale<< ::std::endl;
	::std::cout << "center: " << (sb.min() + sb.max()) / 2.0 << ::std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	void createGround(Geometry::Scene & scene)
///
/// \brief	Adds a sirface area light to the scene
///
/// \author	F. Lamarche, Universit� de Rennes 1
/// \date	03/12/2013
///
/// \param [in,out]	scene	The scene.
////////////////////////////////////////////////////////////////////////////////////////////////////
void createSurfaceLigth(Geometry::Scene & scene, double value)
{
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Geometry::Material * material = new Geometry::Material(RGBColor(), RGBColor(), RGBColor(), 100.0f, RGBColor(value,value,value));
	Geometry::Square square(material);
	Math::Vector3f scaleV = (sb.max() - sb.min());
	//double scale = ::std::max(scaleV[0], scaleV[1])*0.1;
	double factor = 0.5;
	square.scaleX(scaleV[0] * factor);
	square.scaleY(scaleV[1] * factor);
	Math::Vector3f center = (sb.min() + sb.max()) / 2.0;
	center[2] = sb.max()[2]*3;
	square.translate(center);
	scene.add(square);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	void initDiffuse(Geometry::Scene & scene)
///
/// \brief	Adds a Cornell Box with diffuse material on each walls to the scene. This Cornel box
/// 		contains two cubes.
///
/// \author	F. Lamarche, Universit� de Rennes 1
/// \date	03/12/2013
///
/// \param [in,out]	scene	The scene.
////////////////////////////////////////////////////////////////////////////////////////////////////
void initDiffuse(Geometry::Scene & scene)
{
	Geometry::Material * material = new Geometry::Material(RGBColor(), RGBColor(0,0,0.0), RGBColor(0.95f,0.95f,0.95f), 1, RGBColor()) ;
	Geometry::Material * material2 = new Geometry::Material(RGBColor(), RGBColor(1.0,1.0,1.0), RGBColor(0,0,0), 1000, RGBColor()) ;
	Geometry::Material * cubeMat = new Geometry::Material(RGBColor(), RGBColor(1, 0.0, 0.0), RGBColor(0.0, 0.0, 0.0), 20.0f, RGBColor());
	//Geometry::Material * cubeMat = new Geometry::Material(RGBColor(), RGBColor(1.0f,0.0,0.0), RGBColor(0.0,0.0,0.0), 20.0f, RGBColor(10.0,0,0)) ;
	Geometry::Cornel geo(material2, material2, material2, material2, material2, material2) ;

	geo.scaleX(10) ;
	geo.scaleY(10) ;
	geo.scaleZ(10) ;
	scene.add(geo) ;

	Geometry::Cube tmp(cubeMat) ;
	tmp.translate(Math::makeVector(1.5,-1.5,0.0)) ;
	scene.add(tmp) ;

	Geometry::Cube tmp2(cubeMat) ;
	tmp2.translate(Math::makeVector(2,1,-4)) ;
	scene.add(tmp2) ;

	// 2.2 Adds point lights in the scene
	{
		Geometry::PointLight pointLight(Math::makeVector(0.0f, 0.f, 2.0f), RGBColor(0.5f, 0.5f, 0.5f));
		scene.add(pointLight);
	}
	{
		Geometry::PointLight pointLight2(Math::makeVector(4.f, 0.f, 0.f), RGBColor(0.5f, 0.5f, 0.5f));
		scene.add(pointLight2);
	}
	{
		Geometry::Camera camera(Math::makeVector(-4.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f), 0.3f, 1.0f, 1.0f);
		scene.setCamera(camera);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	void initSpecular(Geometry::Scene & scene)
///
/// \brief	Adds a Cornel box in the provided scene. Walls are specular and the box contains two
/// 		cubes.
///
/// \author	F. Lamarche, Universit� de Rennes 1
/// \date	03/12/2013
////////////////////////////////////////////////////////////////////////////////////////////////////
void initSpecular(Geometry::Scene & scene)
{
	Geometry::Material * material = new Geometry::Material(RGBColor(), RGBColor(0,0,0.0), RGBColor(0.7f,0.7f,0.7f), 100, RGBColor()) ;
	Geometry::Material * material2 = new Geometry::Material(RGBColor(), RGBColor(0,0,1.0f), RGBColor(0,0,0), 1000, RGBColor()) ;
	//Geometry::Material * cubeMat = new Geometry::Material(RGBColor(), RGBColor(1.0f,0.0,0.0), RGBColor(0.0,0.0,0.0), 20.0f, RGBColor(10.0,0,0)) ;
	Geometry::Material * cubeMat = new Geometry::Material(RGBColor(), RGBColor(1.0f, 0.0, 0.0), RGBColor(0.0, 0.0, 0.0), 20.0f, RGBColor());
	Geometry::Cornel geo(material, material, material, material, material, material) ; //new Geometry::Cube(material2) ;////new Cone(4, material) ; //new Geometry::Cylinder(5, 1, 1, material) ;////////new Geometry::Cube(material) ;////; //new Geometry::Cube(material) ; //new Geometry::Cylinder(100, 2, 1, material) ; //

	geo.scaleX(10) ;
	geo.scaleY(10) ;
	geo.scaleZ(10) ;
	scene.add(geo) ;

	Geometry::Cube tmp(cubeMat) ;
	tmp.translate(Math::makeVector(1.5,-1.5,0.0)) ;
	scene.add(tmp) ;

	Geometry::Cube tmp2(cubeMat) ;
	tmp2.translate(Math::makeVector(2,1,-4)) ;
	scene.add(tmp2) ;

	// 2.2 Adds point lights in the scene
	{
		Geometry::PointLight pointLight(Math::makeVector(0.0f, 0.f, 2.0f), RGBColor(0.5f, 0.5f, 0.5f)*5);
		scene.add(pointLight);
	}
	{
		Geometry::PointLight pointLight2(Math::makeVector(4.f, 0.f, 0.f), RGBColor(0.5f, 0.5f, 0.5f)*5);
		scene.add(pointLight2);
	}
	// Sets the camera
	{
		Geometry::Camera camera(Math::makeVector(-4.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f), 0.3f, 1.0f, 1.0f);
		scene.setCamera(camera);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	void initDiffuseSpecular(Geometry::Scene & scene)
///
/// \brief	Adds a Cornel box in the provided scene. The cornel box as diffuse and specular walls and
/// 		contains two boxes.
///
/// \author	F. Lamarche, Universit� de Rennes 1
/// \date	03/12/2013
///
/// \param [in,out]	scene	The scene.
////////////////////////////////////////////////////////////////////////////////////////////////////
void initDiffuseSpecular(Geometry::Scene & scene)
{
	Geometry::Material * material = new Geometry::Material(RGBColor(), RGBColor(0,0,0.0), RGBColor(0.7f,0.7f,0.7f), 100, RGBColor()) ;
	Geometry::Material * material2 = new Geometry::Material(RGBColor(), RGBColor(1,1,1.0f), RGBColor(0,0,0), 1000, RGBColor()) ;
	Geometry::Material * cubeMat = new Geometry::Material(RGBColor(), RGBColor(1.0f, 0.0, 0.0), RGBColor(0.0, 0.0, 0.0), 20.0f, RGBColor());
	Geometry::Material * cubeMat2 = new Geometry::Material(RGBColor(), RGBColor(1.0f, 0.0, 0.0), RGBColor(0.0, 0.0, 0.0), 20.0f, RGBColor());
	//Geometry::Material * cubeMat = new Geometry::Material(RGBColor(), RGBColor(0.0f,0.0,0.0), RGBColor(0.0,0.0,0.0), 20.0f, RGBColor(10.0,0,0)) ;
	//Geometry::Material * cubeMat2 = new Geometry::Material(RGBColor(), RGBColor(0.0f,0.0,0.0), RGBColor(0.0,0.0,0.0), 20.0f, RGBColor(0.0,10,0)) ;
	Geometry::Cornel geo(material2, material2, material, material, material, material) ; //new Geometry::Cube(material2) ;////new Cone(4, material) ; //new Geometry::Cylinder(5, 1, 1, material) ;////////new Geometry::Cube(material) ;////; //new Geometry::Cube(material) ; //new Geometry::Cylinder(100, 2, 1, material) ; //

	geo.scaleX(10) ;
	geo.scaleY(10) ;
	geo.scaleZ(10) ;
	scene.add(geo) ;


	Geometry::Cube tmp(cubeMat2) ;
	tmp.translate(Math::makeVector(1.5,-1.5,0.0)) ;
	scene.add(tmp) ;

	Geometry::Cube tmp2(cubeMat) ;
	tmp2.translate(Math::makeVector(2,1,-4)) ;
	scene.add(tmp2) ;

	// 2.2 Adds point lights in the scene
	{
		Geometry::PointLight pointLight(Math::makeVector(0.0f, 0.f, 2.0f), RGBColor(0.5f, 0.5f, 0.5f));
		scene.add(pointLight);
	}
	{
		Geometry::PointLight pointLight2(Math::makeVector(4.f, 0.f, 0.f), RGBColor(0.5f, 0.5f, 0.5f));
		scene.add(pointLight2);
	}
	{
		Geometry::Camera camera(Math::makeVector(-4.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f), 0.3f, 1.0f, 1.0f);
		scene.setCamera(camera);
	}
}

/// <summary>
/// Intializes a scene containing a garage.
/// </summary>
/// <param name="scene"></param>
void initGarage(Geometry::Scene & scene)
{
	Geometry::Loader3ds loader(m_modelDirectory+"\\garage.3ds", m_modelDirectory+"");

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		scene.add(*loader.getMeshes()[cpt]);
	}

	// 2.2 Adds point lights in the scene
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0)*1000);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0)*1000);
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(750.0f, -1500.f, 1000.f)*0.85f, Math::makeVector(200.0f, 0.0f, 0.0f), 0.3f, 1.0f, 1.0f);
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a guitar.
/// </summary>
/// <param name="scene"></param>
void initGuitar(Geometry::Scene & scene)
{
	Geometry::Loader3ds loader(m_modelDirectory+"\\guitar2.3ds", m_modelDirectory+"");

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		scene.add(*loader.getMeshes()[cpt]);
	}
	// 2.2 Adds point lights in the scene
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position + Math::makeVector(0.f, 0.f, 70.f/*100.f*/), RGBColor(1.0, 1.0, 1.0)*500);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position+Math::makeVector(0.f,0.f,200.f), RGBColor(1.0, 1.0, 1.0)*500);
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(-500., -1000., 1000.)*1.05, Math::makeVector(500.f, 0.0f, 0.0f), 0.6f, 1.0f, 1.0f);
		camera.translateLocal(Math::makeVector(100.0f, -100.f, -200.f));
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a dog
/// </summary>
/// <param name="scene"></param>
void initDog(Geometry::Scene & scene)
{
	Geometry::Loader3ds loader(m_modelDirectory+"\\Dog\\dog.3ds", m_modelDirectory+"\\dog");

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		scene.add(*loader.getMeshes()[cpt]);
	}

	// 2.2 Adds point lights in the scene
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0)*2);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0)*2);
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(10.f, 10.f, 6.f)*0.5, Math::makeVector(0.f, 0.0f, 2.5f), .7f, 1.0f, 1.0f);
		camera.translateLocal(Math::makeVector(-0.4f, 0.f, 0.9f));
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a temple.
/// </summary>
/// <param name="scene"></param>
void initTemple(Geometry::Scene & scene)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f,0.0f,0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory+"\\Temple\\Temple of St Seraphim of Sarov N270116_2.3ds", m_modelDirectory+"\\Temple");

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		//loader.getMeshes()[cpt]->translate(Math::makeVector(20.f, 0.f, 40.0f));
		scene.add(*loader.getMeshes()[cpt]);
		box.update(*loader.getMeshes()[cpt]);
	}

	// 2.2 Adds point lights in the scene
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	//Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0)*60.0);
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0)*600.0);
	scene.add(light1);
	position[1] = -position[1];
	//Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0)*30);
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) * 300);
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(20.0f, -100.0f, 15.0f), Math::makeVector(-20.f, 0.f, -40.f), 0.3f, 1.0f, 1.0f);
		camera.translateLocal(Math::makeVector(40.f, 0.f, 0.f));
		scene.setCamera(camera);
	}
	createGround(scene);
	createSurfaceLigth(scene, 50);
	//createSurfaceLigth(scene, 500);
}

/// <summary>
/// Initializes a scene containing a robot.
/// </summary>
/// <param name="scene"></param>
void initRobot(Geometry::Scene & scene)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory+"\\Robot.3ds", m_modelDirectory+"");

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		//loader.getMeshes()[cpt]->translate(Math::makeVector(20.f, 0.f, 40.0f));
		scene.add(*loader.getMeshes()[cpt]);
		box.update(*loader.getMeshes()[cpt]);
	}

	// 2.2 Adds point lights in the scene
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0)*60.0);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) * 30);
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(100.0f, -50.0f, 0.0f), Math::makeVector(0.f, 0.f, -20.f), 0.3f, 1.0f, 1.0f);
		camera.translateLocal(Math::makeVector(0.f, 40.f, 50.f));
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a grave stone
/// </summary>
/// <param name="scene"></param>
void initGraveStone(Geometry::Scene & scene)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory+"\\gravestone\\GraveStone.3ds", m_modelDirectory+"\\gravestone");

	::std::vector<Geometry::Material*> materials = loader.getMaterials();
	for (auto it = materials.begin(), end = materials.end(); it != end; ++it)
	{
		//(*it)->setSpecular(RGBColor());
		//(*it)->setSpecular((*it)->getSpecular()*0.05);
	}

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		//loader.getMeshes()[cpt]->translate(Math::makeVector(20.f, 0.f, 40.0f));
		scene.add(*loader.getMeshes()[cpt]);
		box.update(*loader.getMeshes()[cpt]);
	}

	// 2.2 Adds point lights in the scene
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0)*1500*0.2);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) * 1000*0.4);
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(0.f, -300.0f, 200.0f), Math::makeVector(0.f, 0.f, 60.f), 0.3f, 1.0f, 1.0f);
		camera.translateLocal(Math::makeVector(0.f, 80.f, 120.f));
		scene.setCamera(camera);
	}
	createGround(scene);
	//createSurfaceLigth(scene, 400);
}

/// <summary>
/// Initializes a scene containing a boat
/// </summary>
/// <param name="scene"></param>
void initBoat(Geometry::Scene & scene)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory+"\\Boat\\boat.3ds", m_modelDirectory+"\\boat");
	//// We remove the specular components of the materials...
	::std::vector<Geometry::Material*> materials = loader.getMaterials();
	for (auto it = materials.begin(), end = materials.end(); it != end; ++it)
	{
		(*it)->setSpecular(RGBColor());
	}

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		//loader.getMeshes()[cpt]->translate(Math::makeVector(20.f, 0.f, 40.0f));
		scene.add(*loader.getMeshes()[cpt]);
		box.update(*loader.getMeshes()[cpt]);
	}

	// 2.2 Adds point lights in the scene
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0)*6000.0);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) * 20000);
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(5000.f, 5000.f, 200.0f), Math::makeVector(0.f, 0.f, 60.f), 0.3f, 1.0f, 1.0f);
		camera.translateLocal(Math::makeVector(2000.f, 3000.f, 700.f));
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a tibet house
/// </summary>
/// <param name="scene"></param>
void initTibetHouse(Geometry::Scene & scene)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory+"\\TibetHouse\\TibetHouse.3ds", m_modelDirectory+"\\TibetHouse");
	// We remove the specular components of the materials... and add surface light sources :)
	::std::vector<Geometry::Material*> materials = loader.getMaterials();
	for (auto it = materials.begin(), end = materials.end(); it != end; ++it)
	{
		(*it)->setSpecular(RGBColor());
		if ((*it)->getTextureFile() == m_modelDirectory + "\\TibetHouse" + "\\3D69C2DE.png")
		{
			(*it)->setEmissive(RGBColor(1.0, 1.0, 1.0)*50.0);
		}
	}

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		//loader.getMeshes()[cpt]->translate(Math::makeVector(20.f, 0.f, 40.0f));
		scene.add(*loader.getMeshes()[cpt]);
		box.update(*loader.getMeshes()[cpt]);
	}

	// 2.2 Adds point lights in the scene
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0)*50.0);
	//scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) * 200);
	//scene.add(light2);
	//Geometry::PointLight light3(Math::makeVector(5.f, 35.f, 5.f), RGBColor(1.0, 1.0, 1.0) * 200); //*50
	//scene.add(light3);
	{
		Geometry::Camera camera(Math::makeVector(20.f, 0.f, 0.0f), Math::makeVector(5.f, 35.f, 0.f), 0.3f, 1.0f, 1.0f);
		camera.translateLocal(Math::makeVector(0.f, 5.f, 0.f)/*+Math::makeVector(0.0,5.0,0.0)*/);
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a tibet house. Camera is placed inside the house.
/// </summary>
/// <param name="scene"></param>
void initTibetHouseInside(Geometry::Scene & scene)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory+"\\TibetHouse\\TibetHouse.3ds", m_modelDirectory+"\\TibetHouse");
	// We remove the specular components of the materials... and add surface light sources :)
	::std::vector<Geometry::Material*> materials = loader.getMaterials();
	for (auto it = materials.begin(), end = materials.end(); it != end; ++it)
	{
		(*it)->setSpecular(RGBColor());
		if ((*it)->getTextureFile() == m_modelDirectory + "\\TibetHouse" + "\\3D69C2DE.png")
		{
			(*it)->setEmissive(RGBColor(1.0, 1.0, 1.0)*500.0);
		}
	}

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		//loader.getMeshes()[cpt]->translate(Math::makeVector(20.f, 0.f, 40.0f));
		scene.add(*loader.getMeshes()[cpt]);
		box.update(*loader.getMeshes()[cpt]);
	}

	// 2.2 Adds point lights in the scene
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0)*500.0);
	//scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) * 200);
	//scene.add(light2);
	Geometry::PointLight light3(Math::makeVector(5.f, 35.f, 5.f), RGBColor(1.0, 1.0, 1.0) * 50);
	//scene.add(light3);
	{
		Geometry::Camera camera(Math::makeVector(20.f, 0.f, 5.0f), Math::makeVector(5.f, 35.f, 5.f), 0.3f, 1.0f, 1.0f);
		camera.translateLocal(Math::makeVector(10.f-2, 30.f, 0.f));
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a medieval city
/// </summary>
/// <param name="scene"></param>
void initMedievalCity(Geometry::Scene & scene)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory+"\\Medieval\\MedievalCity.3ds", m_modelDirectory+"\\Medieval\\texture");
	// We remove the specular components of the materials...
	::std::vector<Geometry::Material*> materials = loader.getMaterials();
	for (auto it = materials.begin(), end = materials.end(); it != end; ++it)
	{
		(*it)->setSpecular(RGBColor());
	}

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		//loader.getMeshes()[cpt]->translate(Math::makeVector(20.f, 0.f, 40.0f));
		//loader.getMeshes()[cpt]->rotate(Math::Quaternion<double>(Math::makeVector(0.0, 1.0, 0.0), Math::pi));
		scene.add(*loader.getMeshes()[cpt]);
		box.update(*loader.getMeshes()[cpt]);
	}

	// 2.2 Adds point lights in the scene
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position+Math::makeVector(0.0,0.0,150.0), RGBColor(1.0, 0.6, 0.3)*800.0);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position + Math::makeVector(0.0, 0.0, 1000.0), RGBColor(1.0, 1.0, 1.0) * 400);
	scene.add(light2);
	//Geometry::PointLight light3(Math::makeVector(5.f, 35.f, 5.f), RGBColor(1.0, 1.0, 1.0) * 50);
	//scene.add(light3);
	{
		Geometry::Camera camera(Math::makeVector(0.f, 300.f, 1000.0f), Math::makeVector(0.0,0.0,0.0), 0.3f, 1.0f, 1.0f);
		camera.translateLocal(Math::makeVector(0.0, 800., -100.0));
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a sombrero
/// </summary>
/// <param name="scene"></param>
void initSombrero(Geometry::Scene & scene)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory+"/sombrero/sombrero.3ds", m_modelDirectory+"/sombrero");
	// We remove the specular components of the materials...
	::std::vector<Geometry::Material*> materials = loader.getMaterials();
	for (auto it = materials.begin(), end = materials.end(); it != end; ++it)
	{
		(*it)->setSpecular(RGBColor());
	}

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		scene.add(*loader.getMeshes()[cpt]);
		box.update(*loader.getMeshes()[cpt]);
	}

	// 2.2 Adds point lights in the scene
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0)*50.0);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) * 50);
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(300.f, 0.f, 100.0f), Math::makeVector(0.f, 0.f, 0.f), 0.3f, 1.0f, 1.0f);
		//camera.translateLocal(Math::makeVector(2000.f, 3000.f, 700.f));
		scene.setCamera(camera);
	}
	createGround(scene);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	void waitKeyPressed()
///
/// \brief	Waits until a key is pressed.
///
/// \author	F. Lamarche, Universit� de Rennes 1
/// \date	03/12/2013
////////////////////////////////////////////////////////////////////////////////////////////////////
void waitKeyPressed()
{
  SDL_Event event;
  bool done = false;
  while( !done ) {
    while ( SDL_PollEvent(&event) ) {
      switch (event.type) {
        case SDL_KEYDOWN:
        /*break;*/
        case SDL_QUIT:
          done = true;
        break;
        default:
        break;
      }
    }/*while*/
  }/*while(!done)*/
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	int main(int argc, char ** argv)
///
/// \brief	Main entry-point for this application.
///
/// \author	F. Lamarche, Universit� de Rennes 1
/// \date	03/12/2013
///
/// \param	argc	Number of command-line arguments.
/// \param	argv	Array of command-line argument strings.
///
/// \return	Exit-code for the process - 0 for success, else an error code.
////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char ** argv)
{
	omp_set_num_threads(8);

	// 1 - Initializes a window for rendering
	//Visualizer::Visualizer visu(1000,1000) ;
	Visualizer::Visualizer visu(500, 500);
	// Visualizer::Visualizer visu(300,300) ;

	// 2 - Initializes the scene
	Geometry::Scene scene(&visu) ;

	// 2.1 initializes the geometry (choose only one initialization)
	// initDiffuse(scene) ;
	// initDiffuseSpecular(scene) ;
	// initSpecular(scene) ;
	//initGuitar(scene);
	//initDog(scene);
	//initGarage(scene);
	//initTemple(scene);
	//initRobot(scene);
	//initGraveStone(scene);
	//initBoat(scene);
	initSombrero(scene);
	//initTibetHouse(scene);
	//initTibetHouseInside(scene);
	//initMedievalCity(scene);

	// Shows stats
	scene.printStats();

	// 3 - Computes the scene
	// unsigned int passPerPixel = 64/8;	// Number of rays per pixel
	unsigned int passPerPixel = 1;	// Number of rays per pixel

	unsigned int subPixelSampling = 4;	// Antialiasing
	//unsigned int maxBounce = 5;
	unsigned int maxBounce = 5;			// Maximum number of bounces

	//scene.setDiffuseSamples(16);
	//scene.setSpecularSamples(16);
	scene.setDiffuseSamples(1);
	scene.setSpecularSamples(1);
	//scene.setDiffuseSamples(32);
	//scene.setSpecularSamples(32);
	//scene.setDiffuseSamples(16);
	//scene.setSpecularSamples(16);
	//scene.setDiffuseSamples(4);
	//scene.setSpecularSamples(4);
	scene.setLightSamples(32);
	scene.compute(maxBounce, subPixelSampling, passPerPixel) ;

	// 4 - waits until a key is pressed
	waitKeyPressed() ;

	return 0 ;
}
