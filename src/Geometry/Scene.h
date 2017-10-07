#ifndef _Geometry_Scene_H
#define _Geometry_Scene_H

//#include <windows.h>
#include <Geometry/Geometry.h>
#include <Geometry/PointLight.h>
#include <Visualizer/Visualizer.h>
#include <Geometry/Camera.h>
#include <Geometry/BoundingBox.h>
#include <Math/RandomDirection.h>
//#include <windows.h>
#include <System/aligned_allocator.h>
#include <Math/Constant.h>
#include <queue>
#include <functional>
#include <random>
#include <Geometry/LightSampler.h>
#include <math.h>

namespace Geometry
{
	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \class	Scene
	///
	/// \brief	An instance of a geometric scene that can be rendered using ray casting. A set of methods
	/// 		allowing to add geometry, lights and a camera are provided. Scene rendering is achieved by
	/// 		calling the Scene::compute method.
	///
	/// \author	F. Lamarche, Universit� de Rennes 1
	/// \date	03/12/2013
	////////////////////////////////////////////////////////////////////////////////////////////////////
	class Scene
	{
	protected:
		/// \brief	The visualizer (rendering target).
		Visualizer::Visualizer * m_visu ;
		/// \brief	The scene geometry (basic representation without any optimization).
		::std::deque<::std::pair<BoundingBox, Geometry> > m_geometries ;
		//Geometry m_geometry ;
		/// \brief	The lights.
		std::vector<PointLight> m_lights ;
		/// \brief	The camera.
		Camera m_camera ;
		/// \brief The scene bounding box
		BoundingBox m_sceneBoundingBox;
		/// \brief number of diffuse samples for global illumination
		size_t m_diffuseSamples ;
		/// \brief Number of specular samples
		size_t m_specularSamples ;
		/// \brief Number of light samples
		size_t m_lightSamples;
		/// brief Rendering pass number
		int m_pass;
		/// <summary>
		/// The light sampler associated with the scene
		/// </summary>
		LightSampler m_lightSampler;



	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Scene::Scene(Visualizer::Visualizer * visu)
		///
		/// \brief	Constructor.
		///
		/// \author	F. Lamarche, Universit� de Rennes 1
		/// \date	03/12/2013
		///
		/// \param [in,out]	visu	If non-null, the visu.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Scene(Visualizer::Visualizer * visu)
			: m_visu(visu), m_diffuseSamples(30), m_specularSamples(30), m_lightSamples(0)
		{}

		/// <summary>
		/// Prints stats about the geometry associated with the scene
		/// </summary>
		void printStats()
		{
			size_t nbTriangles = 0;
			for (auto it = m_geometries.begin(), end = m_geometries.end(); it != end; ++it)
			{
				nbTriangles += it->second.getTriangles().size();
			}
			::std::cout << "Scene: " << nbTriangles << " triangles" << ::std::endl;
		}

		/// <summary>
		/// Computes the scene bounding box.
		/// </summary>
		/// <returns></returns>
		const BoundingBox & getBoundingBox()
		{
			return m_sceneBoundingBox;
		}

		/// <summary>
		/// Sets the number of diffuse samples
		/// </summary>
		/// <param name="number"> The number of diffuse samples</param>
		void setDiffuseSamples(size_t number)
		{
			m_diffuseSamples = number;
		}

		/// <summary>
		/// Sets the number of specular samples
		/// </summary>
		/// <param name="number"></param>
		void setSpecularSamples(size_t number)
		{
			m_specularSamples = number;
		}

		/// <summary>
		/// Sets the number of light samples if the scene contains surface area lights
		/// </summary>
		/// <param name="number">The number of samples</param>
		void setLightSamples(size_t number)
		{
			m_lightSamples = number;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Scene::add(const Geometry & geometry)
		///
		/// \brief	Adds a geometry to the scene.
		///
		/// \author	F. Lamarche, Universit� de Rennes 1
		/// \date	03/12/2013
		///
		/// \param	geometry The geometry to add.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void add(const Geometry & geometry)
		{
			if (geometry.getVertices().size() == 0) { return; }
			BoundingBox box(geometry) ;
			m_geometries.push_back(::std::make_pair(box, geometry)) ;
			m_geometries.back().second.computeVertexNormals(Math::piDiv4/2);
			if (m_geometries.size() == 1)
			{
				m_sceneBoundingBox = box;
			}
			else
			{
				m_sceneBoundingBox.update(box);
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Scene::add(PointLight * light)
		///
		/// \brief	Adds a poitn light in the scene.
		///
		/// \author	F. Lamarche, Universit� de Rennes 1
		/// \date	04/12/2013
		///
		/// \param [in,out]	light	If non-null, the light to add.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void add(const PointLight & light)
		{
			m_lights.push_back(light) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Scene::setCamera(Camera const & cam)
		///
		/// \brief	Sets the camera.
		///
		/// \author	F. Lamarche, Universit� de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	cam	The camera.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void setCamera(Camera const & cam)
		{
			m_camera = cam ;
		}

		RayTriangleIntersection findIntersection(CastedRay const & cRay){
			//Methode naive, parcours de tous les triangles de la sc�ne
			for(auto it = this->m_geometries.begin(); it != this->m_geometries.end(); ++it){
				//parcours de toutes les g�om�tries
				//memorisation de  la geo courante
				Geometry & geoTmp((*it).second);

				//parcours de tous les triangles de la g�om�trie
				for(auto itTriangle = geoTmp.getTriangles().begin(); itTriangle != geoTmp.getTriangles().end(); ++itTriangle){
						//on calcule l'intersection (la + proche est m�moris�)
						cRay.intersect(&(*itTriangle));
				}
			}

			if(cRay.validIntersectionFound()){ //if we found an intersection, return it
				return cRay.intersectionFound();
			}
			else{ //otherwise return an unvalid intersection
				return RayTriangleIntersection();
			}
		}

		bool inShadow(PointLight const & pl, RayTriangleIntersection const & rTI){
			bool shadow = false;
			Math::Vector3f lightToIntersection(rTI.intersection()-pl.position());
			CastedRay shadowRay(pl.position(),lightToIntersection);

			RayTriangleIntersection shadowRayIntersection = findIntersection(shadowRay);
			if(shadowRayIntersection.valid()){
				if(shadowRayIntersection.triangle()!=rTI.triangle()){ //on regarde si on intersecte le même triangle
					shadow = true; //l'objet n'est pas éclairé
				}
			}
			return shadow;
		}

		Math::Vector3f getReflectionRay(PointLight const & pl, RayTriangleIntersection const & rTI, Math::Vector3f const & cameraOrientedNormale /*name it N*/){
			//first step, the reflected ray
			Math::Vector3f intersectToLight(pl.position()-rTI.intersection()); //name it L
			//the vector should be normalided to compute the cos. Camera is already normalized
			 //the formule is (2(N.L)).N-L
			Math::Vector3f reflexion((cameraOrientedNormale*(2.0*(cameraOrientedNormale*intersectToLight.normalized())))-intersectToLight.normalized());
			return reflexion;
		}

		RGBColor computeSpecularColor(PointLight const & pl, RayTriangleIntersection const & rTI, Math::Vector3f const & cameraOrientedNormale, Math::Vector3f const & intersectToCamera){
			//the first step is to compute the reflected ray. We name it R
			Math::Vector3f reflexion(getReflectionRay(pl,rTI,cameraOrientedNormale));

			/*
			* Lets consider
			* - K_s for the specular triangle color
			* - I_l for the light color
			* - V for the vector which goes from intersection to camera
			* - n_shininess for the shininess coef of the triangle
			* The formula is: I_spec = K_s * I_l (V.R)^{n_shininess}
			*/
			RGBColor specular(rTI.triangle()->material()->getSpecular()*pl->color())*pow(reflexion*intersectToCamera.normalized(),rTI.triangle()->material()->getShininess());

			return specular;
		}

		RGBColor computeColor(RayTriangleIntersection const & rTI, CastedRay const & ray){
			Triangle * triangle = rTI.triangle();
			Math::Vector3f lightDirection(0.0);
			Math::Vector3f lightPosition(0.0);
			Math::Vector3f invertLightDirection(0.0);

			Math::Vector3f intersection = rTI.intersection();
			Math::Vector3f normale(triangle->normal(intersection));
			RGBColor colorResult(0.0, 0.0, 0.0);
			RGBColor colorSpecular(0.0,0.0,0.0);
			RGBColor colorDiffuse(0.0,0.0,0.0);


			float scalarNL = 0;

			//on construit le vecteur qui va de l'intersection vers la source
			Math::Vector3f intersectToSource(ray.source()-intersection);

			//on regarde si la normale est dans le bon sens (produit scalaire positif)
			if(normale*intersectToSource<0){
				normale = Math::Vector3f(0.0)-normale;
			}

			//on parcourt toutes les sources de lumière et on regarde la couleur:
			for(auto itLight = this->m_lights.begin(); itLight != this->m_lights.end(); ++itLight){
				lightPosition = itLight->position(); //on prends la position de la lumière (points)
				//on part de l'arrivée vers le départ (lampe vers l'intersection), construction du vecteur
				lightDirection = lightPosition-intersection;

				//on regarde si la source vient de devant ou derrière (comparaison avec la normale produit scalaire)
				scalarNL = lightDirection.normalized()*normale.normalized(); //le scalaire à la normale, pour être comparable avec cette dernière doit être nomalisé
				if(scalarNL>0){ // si la source est dans notre demi espace (elle éclaire alors les objets devant nous)
					if(!inShadow(*itLight,rTI)){ //est ce que la lumière éclaire bien notre objet?
						colorDiffuse = (triangle->material()->getDiffuse()*itLight->color())*scalarNL;
						colorSpecular = computeSpecularColor(*itLight,rTI,normale,intersectToSource);
						colorResult = colorResult + colorDiffuse + colorSpecular;
					}
				}
			}
			return colorResult;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	RGBColor Scene::sendRay(Ray const & ray, double limit, int depth, int maxDepth)
		///
		/// \brief	Sends a ray in the scene and returns the computed color
		///
		/// \author	F. Lamarche, Universit� de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	ray			The ray.
		/// \param	depth   	The current depth.
		/// \param	maxDepth	The maximum depth.
		///
		/// \return	The computed color.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		RGBColor sendRay(Ray const & ray, int depth, int maxDepth, int diffuseSamples, int specularSamples)
		{
			RGBColor result(0.0, 0.0, 0.0);
			CastedRay cRay(ray);

			RayTriangleIntersection primaryIntersection = findIntersection(cRay);

			if(primaryIntersection.valid()){ //if we have a triangle
				//compute the color with light
				result = computeColor(primaryIntersection,cRay);
			}

			return result;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Scene::compute(int maxDepth)
		///
		/// \brief	Computes a rendering of the current scene, viewed by the camera.
		///
		/// \author	F. Lamarche, Universit� de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	maxDepth	The maximum recursive depth.
		/// \param  subPixelDivision subpixel subdivisions to handle antialiasing
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void compute(int maxDepth, int subPixelDivision = 1, int passPerPixel = 1)
		{
			// We prepare the light sampler (the sampler only stores triangles with a non null emissive component).
			for (auto it = m_geometries.begin(), end = m_geometries.end(); it != end; ++it)
			{
				m_lightSampler.add(it->second);
			}

			// Step on x and y for subpixel sampling
			double step = 1.0f/subPixelDivision ;
			// Table accumulating values computed per pixel (enable rendering of each pass)
			::std::vector<::std::vector<::std::pair<int, RGBColor> > > pixelTable(m_visu->width(), ::std::vector<::std::pair<int, RGBColor> >(m_visu->width(), ::std::make_pair(0, RGBColor()))) ;

			// 1 - Rendering time
			// LARGE_INTEGER frequency;        // ticks per second
			// LARGE_INTEGER t1, t2;           // ticks
			double elapsedTime;
			// get ticks per second
			// QueryPerformanceFrequency(&frequency);
			// start timer
			// QueryPerformanceCounter(&t1);
			// Rendering pass number
			m_pass = 0;
			// Rendering
			for(int passPerPixelCounter = 0 ; passPerPixelCounter<passPerPixel ; ++passPerPixelCounter)
			{
				for (double xp = -0.5; xp < 0.5; xp += step)
				{
					for (double yp = -0.5; yp < 0.5; yp += step)
					{
						::std::cout << "Pass: " << m_pass << "/" << passPerPixel * subPixelDivision * subPixelDivision << ::std::endl;
						++m_pass;
						// Sends primary rays for each pixel (uncomment the pragma to parallelize rendering)
#pragma omp parallel for schedule(dynamic)//, 10)//guided)//dynamic)
						for (int y = 0; y < m_visu->height(); y++)
						{
							for (int x = 0; x < m_visu->width(); x++)
							{
#pragma omp critical (visu)
								m_visu->plot(x, y, RGBColor(1000.0, 0.0, 0.0));
								// Ray casting
								RGBColor result = sendRay(m_camera.getRay(((double)x + xp) / m_visu->width(), ((double)y + yp) / m_visu->height()), 0, maxDepth, m_diffuseSamples, m_specularSamples);
								// Accumulation of ray casting result in the associated pixel
								::std::pair<int, RGBColor> & currentPixel = pixelTable[x][y];
								currentPixel.first++;
								currentPixel.second = currentPixel.second + result;
								// Pixel rendering (with simple tone mapping)
#pragma omp critical (visu)
								m_visu->plot(x, y, pixelTable[x][y].second / (double)(pixelTable[x][y].first) * 10);
								// Updates the rendering context (per pixel) - warning per pixel update can be costly...
//#pragma omp critical (visu)
								//m_visu->update();
							}
							// Updates the rendering context (per line)
#pragma omp critical (visu)
							m_visu->update();
						}
						// Updates the rendering context (per pass)
						//m_visu->update();
						// We print time for each pass
						// QueryPerformanceCounter(&t2);
						// elapsedTime = (double)(t2.QuadPart - t1.QuadPart) / (double)frequency.QuadPart;
						double remainingTime = (elapsedTime / m_pass)*(passPerPixel * subPixelDivision * subPixelDivision - m_pass);
						::std::cout << "time: " << elapsedTime << "s. " <<", remaining time: "<< remainingTime << "s. " <<", total time: "<< elapsedTime + remainingTime << ::std::endl;
					}
				}
			}
			// stop timer
			// QueryPerformanceCounter(&t2);
			// elapsedTime = (double)(t2.QuadPart - t1.QuadPart) / (double)frequency.QuadPart;
			::std::cout<<"time: "<<elapsedTime<<"s. "<<::std::endl ;
		}
	} ;
}

#endif
