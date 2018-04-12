#include <DsosUtility.h>

using namespace System::Diagnostics;
using namespace System::IO;

namespace Dsos
{
	bool DsosUtility::CreateIdfFile(OpenStudio::Model^ osModel, String^ idfPathName)
	{
		OpenStudio::EnergyPlusForwardTranslator^ osForwardTranslator = gcnew OpenStudio::EnergyPlusForwardTranslator();
		OpenStudio::Workspace^ osWorkspace = osForwardTranslator->translateModel(osModel);
		OpenStudio::IdfFile^ osIdfFile = osWorkspace->toIdfFile();
		OpenStudio::Path^ osIdfPath = gcnew OpenStudio::Path(idfPathName);
		bool idfSaveCondition = osIdfFile->save(osIdfPath);
		return idfSaveCondition;
	}

	OpenStudio::Model^ DsosUtility::BuildOsModel(
		String^ osmTemplatePath, 
		String^ epwWeatherPath, 
		String^ ddyPath, 
		String^ osmOutputPath, 
		CellComplex^ buildingCellComplex,
		String^ buildingName, 
		String^ buildingType, 
		String^ spaceType, 
		double heatingTemp,
		double coolingTemp,
		double buildingHeight, 
		List<double>^ floorLevels, 
		int numFloors, 
		double glazingRatio, 
		List<Cell^>^ contextBuildings)
	{
		List<OpenStudio::Space^>^ osSpaces = gcnew List<OpenStudio::Space^>();
		OpenStudio::Model^ osModel = GetModelFromTemplate(osmTemplatePath, epwWeatherPath, ddyPath);
		OpenStudio::Building^ osBuilding = ComputeBuilding(osModel, buildingName, buildingType, buildingHeight, numFloors, spaceType);
		List<Cell^>^ pBuildingCells = buildingCellComplex->Cells();
		for each(Cell^ buildingCell in pBuildingCells)
		{
			OpenStudio::Space^ osSpace = AddSpace(
				buildingCell,
				osModel,
				Autodesk::DesignScript::Geometry::Vector::ZAxis(),
				buildingHeight,
				floorLevels,
				glazingRatio,
				heatingTemp,
				coolingTemp
			);

			for each(OpenStudio::Space^ osExistingSpace in osSpaces)
			{
				osSpace->matchSurfaces(osExistingSpace);
			}

			osSpaces->Add(osSpace);
		}

		if (contextBuildings != nullptr)
		{
			for each(Cell^ contextBuilding in contextBuildings)
			{
				AddShadingSurfaces(contextBuilding, osModel);
			}
		}

		osModel->purgeUnusedResourceObjects();

		bool saveCondition = SaveModel(osModel, osmOutputPath);

		if (saveCondition)
		{
			return osModel;
		}

		return nullptr;
	}

	void DsosUtility::PerformEnergyAnalysis(String^ strOsmPath, String^ epwPathName, String^ oswPathName, String^ openStudioExePath)
	{
		OpenStudio::WorkflowJSON^ workflow = gcnew OpenStudio::WorkflowJSON();
		OpenStudio::Path^ osmPath = gcnew OpenStudio::Path(strOsmPath);
		OpenStudio::Path^ osmWeather = gcnew OpenStudio::Path(epwPathName);
		OpenStudio::Path^ osmOswPath = gcnew OpenStudio::Path(oswPathName);
		workflow->setSeedFile(osmPath);
		workflow->setWeatherFile(osmWeather);
		workflow->saveAs(osmOswPath);

		// https://stackoverflow.com/questions/5168612/launch-program-with-parameters
		String^ args = "run -w \"" + oswPathName + "\"";
		System::Diagnostics::ProcessStartInfo^ startInfo = gcnew ProcessStartInfo(openStudioExePath, args);
		startInfo->WorkingDirectory = Path::GetDirectoryName(oswPathName);

		Process^ process = Process::Start(startInfo);
	}

	bool DsosUtility::SaveModel(OpenStudio::Model^ osModel, String^ osmPathName)
	{
		// Purge unused resources
		osModel->purgeUnusedResourceObjects();

		// Create a path string
		OpenStudio::Path^ osPath = gcnew OpenStudio::Path(osmPathName);
		bool osCondition = osModel->save(osPath, true);
		return osCondition;
	}

	OpenStudio::Model^ DsosUtility::GetModelFromTemplate(String^ osmTemplatePath, String^ epwWeatherPath, String^ ddyPath)
	{
		OpenStudio::Path^ osTemplatePath = gcnew OpenStudio::Path(osmTemplatePath);

		// Create an abstract model
		OpenStudio::OptionalModel^ osOptionalModel = OpenStudio::Model::load(osTemplatePath);
		OpenStudio::Model^ osModel = osOptionalModel->__ref__();

		// Read an EPW weather file
		OpenStudio::Path^ osEPWPath = gcnew OpenStudio::Path(epwWeatherPath);
		OpenStudio::EpwFile^ osEPWFile = gcnew OpenStudio::EpwFile(osEPWPath);
		OpenStudio::WeatherFile^ osWeatherFile = osModel->getWeatherFile();
		OpenStudio::WeatherFile::setWeatherFile(osModel, osEPWFile);

		// Read an DDY design days files
		OpenStudio::Path^ osDDYPath = gcnew OpenStudio::Path(ddyPath);
		OpenStudio::EnergyPlusReverseTranslator^ osTranslator = gcnew OpenStudio::EnergyPlusReverseTranslator();
		OpenStudio::OptionalModel^ tempModel01 = osTranslator->loadModel(osDDYPath);
		OpenStudio::Model^ tempModel02 = tempModel01->__ref__();
		OpenStudio::DesignDayVector^ designDays = tempModel02->getDesignDays();
		OpenStudio::DesignDayVector::DesignDayVectorEnumerator^ designDaysEnumerator = designDays->GetEnumerator();

		while (designDaysEnumerator->MoveNext())
		{
			OpenStudio::DesignDay^ aDesignDay = designDaysEnumerator->Current;
			OpenStudio::IdfObject^ anIdfObject = aDesignDay->idfObject();
			osModel->addObject(anIdfObject);
		}
		return osModel;
	}

	OpenStudio::ThermalZone^ DsosUtility::CreateThermalZone(OpenStudio::Model^ model, OpenStudio::Space^ space, double ceilingHeight, double heatingTemp, double coolingTemp)
	{
		//dsosUtility dsos = new dsosUtility();
        
		// Create a thermal zone for the space
		OpenStudio::ThermalZone^ osThermalZone = gcnew OpenStudio::ThermalZone(model);
		osThermalZone->setName(space->name() + "_TZ");
		osThermalZone->setUseIdealAirLoads(true);
		osThermalZone->setCeilingHeight(ceilingHeight);
		osThermalZone->setVolume(space->volume());

		// Assign Thermal Zone to space
		// aSpace.setThermalZone(osThermalZone);//Not available in C#
		OpenStudio::UUID^ tzHandle = osThermalZone->handle();
		int location = 10;
		space->setPointer(location, tzHandle);

		OpenStudio::ScheduleConstant^ heatingScheduleConstant = gcnew OpenStudio::ScheduleConstant(model);
		heatingScheduleConstant->setValue(heatingTemp);
		OpenStudio::ScheduleConstant^ coolingScheduleConstant = gcnew OpenStudio::ScheduleConstant(model);
		coolingScheduleConstant->setValue(coolingTemp);

		// Create a Thermostat
		OpenStudio::ThermostatSetpointDualSetpoint^ osThermostat = gcnew OpenStudio::ThermostatSetpointDualSetpoint(model);

		// Set Heating and Cooling Schedules on the Thermostat
		osThermostat->setHeatingSetpointTemperatureSchedule(heatingScheduleConstant);
		osThermostat->setCoolingSetpointTemperatureSchedule(coolingScheduleConstant);

		// Assign Thermostat to the Thermal Zone
		osThermalZone->setThermostatSetpointDualSetpoint(osThermostat);
		return osThermalZone;
	}

	OpenStudio::BuildingStory^ DsosUtility::AddBuildingStory(OpenStudio::Model^ model, int floorNumber)
	{
		OpenStudio::BuildingStory^ osBuildingStory = gcnew OpenStudio::BuildingStory(model);
		osBuildingStory->setName("STORY_" + floorNumber);
		osBuildingStory->setDefaultConstructionSet(getDefaultConstructionSet(model));
		osBuildingStory->setDefaultScheduleSet(getDefaultScheduleSet(model));
		return osBuildingStory;
	}

	OpenStudio::Building^ DsosUtility::ComputeBuilding(
		OpenStudio::Model^ osModel,
		String^ buildingName,
		String^ buildingType,
		double buildingHeight,
		int numFloors,
		String^ spaceType)
	{
		OpenStudio::Building^ osBuilding = osModel->getBuilding();
		//building.setNumberOfStories(stories);
		osBuilding->setStandardsNumberOfStories(numFloors);
		osBuilding->setDefaultConstructionSet(getDefaultConstructionSet(osModel));
		osBuilding->setDefaultScheduleSet(getDefaultScheduleSet(osModel));
		osBuilding->setName(buildingName);
		osBuilding->setStandardsBuildingType(buildingType);
		double floorToFloorHeight = (double)buildingHeight / (double)numFloors;
		osBuilding->setNominalFloortoFloorHeight(floorToFloorHeight);
		// Get all space types and find the one that matches
		OpenStudio::SpaceTypeVector^ spaceTypes = osModel->getSpaceTypes();
		OpenStudio::SpaceTypeVector::SpaceTypeVectorEnumerator^ spaceTypesEnumerator = spaceTypes->GetEnumerator();
		while (spaceTypesEnumerator->MoveNext())
		{
			OpenStudio::SpaceType^ aSpaceType = spaceTypesEnumerator->Current;
			String^ spaceTypeName = aSpaceType->name()->__str__();
			if (spaceTypeName == spaceType)
			{
				osBuilding->setSpaceType(aSpaceType);
			}
		}
		buildingStories = CreateBuildingStories(osModel, numFloors);
		return osBuilding;
	}

	List<OpenStudio::BuildingStory^>^ DsosUtility::CreateBuildingStories(OpenStudio::Model^ osModel, int numFloors)
	{
		List<OpenStudio::BuildingStory^>^ osBuildingStories = gcnew List<OpenStudio::BuildingStory^>();
		for (int i = 0; i < numFloors; i++)
		{
			osBuildingStories->Add(AddBuildingStory(osModel, (i + 1)));

		}
		return osBuildingStories;
	}

	OpenStudio::DefaultScheduleSet^ DsosUtility::getDefaultScheduleSet(OpenStudio::Model^ model)
	{
		// Get list of default schedule sets
		OpenStudio::DefaultScheduleSetVector^ defaultScheduleSets = model->getDefaultScheduleSets();
		OpenStudio::DefaultScheduleSetVector::DefaultScheduleSetVectorEnumerator^ defSchedEnum = defaultScheduleSets->GetEnumerator();
		defSchedEnum->MoveNext();
		defaultScheduleSet = defSchedEnum->Current;
		return defaultScheduleSet;
	}

	OpenStudio::DefaultConstructionSet^ DsosUtility::getDefaultConstructionSet(OpenStudio::Model ^ model)
	{
		// Get list of default construction sets
		OpenStudio::DefaultConstructionSetVector^ defaultConstructionSets = model->getDefaultConstructionSets();
		// Get the first item and use as the default construction set
		OpenStudio::DefaultConstructionSetVector::DefaultConstructionSetVectorEnumerator^ defConEnum = defaultConstructionSets->GetEnumerator();
		defConEnum->MoveNext();
		defaultConstructionSet = defConEnum->Current;
		return defaultConstructionSet;
	}

	Face^ DsosUtility::ApplyAperture(Face^ face, Face^ apertureDesign, int numEdgeSamples)
	{
		if (numEdgeSamples <= 0)
		{
			throw gcnew Exception("numEdgeSamples must be positive.");
		}
		// 1. Convert the apertures and boundary as faces.
		Wire^ pOuterApertureWire = apertureDesign->OuterBoundary();
		List<Wire^>^ pApertureWires = apertureDesign->InnerBoundaries();

		List<Face^>^ pFaces = gcnew List<Face^>();

		// 2. For each wires, iterate through the edges, sample points, and map them to the 
		for each(Wire^ pApertureWire in pApertureWires)
		{
			List<Edge^>^ pApertureEdges = pApertureWire->Edges();
			List<Edge^>^ pMappedApertureEdges = gcnew List<Edge^>();

			for each(Edge^ pApertureEdge in pApertureEdges)
			{
				List<Vertex^>^ pMappedSampleVertices = gcnew List<Vertex^>();
				for (int i = 0; i < numEdgeSamples; ++i)
				{
					double t = (double)i / (double)(numEdgeSamples-1);
					if (t < 0.0)
					{
						t = 0.0;
					}
					else if (t > 1.0)
					{
						t = 1.0;
					}

					// Find the actual point on the edge
					Vertex^ pSampleVertex = pApertureEdge->PointAtParameter(t);

					// Find the UV-coordinate of the point on the aperture design
					Autodesk::DesignScript::Geometry::UV^ pUV = apertureDesign->UVParameterAtPoint(pSampleVertex);
					double checkedU = pUV->U, checkedV = pUV->V;
					if (checkedU < 0.0)
					{
						checkedU = 0.0;
					}else if (checkedU > 1.0)
					{
						checkedU = 1.0;
					}

					if (checkedV < 0.0)
					{
						checkedV = 0.0;
					}
					else if (checkedV > 1.0)
					{
						checkedV = 1.0;
					}

					pUV = Autodesk::DesignScript::Geometry::UV::ByCoordinates(checkedU, checkedV);

					// Find the point with the same UV-coordinate on the surface, add it to the list
					Vertex^ pMappedSampleVertex = face->PointAtParameter(pUV);
					pMappedSampleVertices->Add(pMappedSampleVertex);
				}

				// Interpolate the mapped vertices to an edge.
				Edge^ pMappedApertureEdge = Edge::ByVertices(pMappedSampleVertices);
				pMappedApertureEdges->Add(pMappedApertureEdge);
			}

			// Connect the mapped edges to a wire
			Wire^ pMappedApertureWire = Wire::ByEdges(pMappedApertureEdges);

			//// Use the wire to make a face on the same supporting surface as the input face's
			Face^ pMappedApertureFace = face->Trim(pMappedApertureWire);
			pFaces->Add(pMappedApertureFace);

			// and attach it as an aperture to the face.
			Context^ pFaceContext = Context::ByTopologyParameters(face, 0.0, 0.0, 0.0);
			Aperture^ pMappedAperture = Aperture::ByTopologyContext(pMappedApertureFace, pFaceContext);
		}
		
		// TODO: should return a copy
		return face;
	}

	OpenStudio::Space^ DsosUtility::AddSpace(
		Cell^ cell, 
		OpenStudio::Model^ osModel, 
		Autodesk::DesignScript::Geometry::Vector^ upVector, 
		double buildingHeight, 
		List<double>^ floorLevels, 
		double glazingRatio, 
		double heatingTemp, 
		double coolingTemp)
	{
		OpenStudio::Space^ osSpace = gcnew OpenStudio::Space(osModel);
		osSpace->setName("Storey"); // todo

		List<Face^>^ faces = cell->Faces();
		List<OpenStudio::Point3dVector^>^ facePointsList = gcnew List<OpenStudio::Point3dVector^>();
		for each(Face^ face in faces)
		{
			OpenStudio::Point3dVector^ facePoints = GetFacePoints(face);
			facePointsList->Add(facePoints);
		}

		for (int i = 0; i < faces->Count; ++i)
		{
			AddSurface(faces[i], facePointsList[i], osSpace, osModel, upVector, glazingRatio);
		}

		int storeyNumber = StoreyNumber(cell, buildingHeight, floorLevels);
		OpenStudio::BuildingStory^ buildingStorey = buildingStories[storeyNumber];
		osSpace->setBuildingStory(buildingStorey);
		osSpace->setDefaultConstructionSet(getDefaultConstructionSet(osModel));
		osSpace->setDefaultScheduleSet(getDefaultScheduleSet(osModel));

		// Get all space types
		OpenStudio::SpaceTypeVector^ osSpaceTypes = osModel->getSpaceTypes();
		OpenStudio::SpaceTypeVector::SpaceTypeVectorEnumerator^ osSpaceTypesEnumerator = osSpaceTypes->GetEnumerator();
		int spaceTypeCount = osSpaceTypes->Count;
		while (osSpaceTypesEnumerator->MoveNext())
		{
			OpenStudio::SpaceType^ osSpaceType = osSpaceTypesEnumerator->Current;
			OpenStudio::OptionalString^ osSpaceTypeOptionalString = osSpaceType->name();
			String^ spaceTypeName = osSpaceTypeOptionalString->__str__();
			if (spaceTypeName == "ASHRAE 189::1-2009 ClimateZone 4-8 MediumOffice")
			{
				osSpace->setSpaceType(osSpaceType);
			}
		}

		Autodesk::DesignScript::Geometry::Solid^ solid = safe_cast<Autodesk::DesignScript::Geometry::Solid^>(cell->Geometry);
		List<Autodesk::DesignScript::Geometry::Geometry^>^ solids = gcnew List<Autodesk::DesignScript::Geometry::Geometry^>();
		solids->Add(solid);
		Autodesk::DesignScript::Geometry::BoundingBox^ boundingBox =
			Autodesk::DesignScript::Geometry::BoundingBox::ByGeometry(solids);
		double ceilingHeight = Math::Abs(boundingBox->MaxPoint->Z - boundingBox->MinPoint->Z);
		OpenStudio::ThermalZone^ thermalZone = CreateThermalZone(osModel, osSpace, ceilingHeight, heatingTemp, coolingTemp);

		return osSpace;
	}


	void DsosUtility::AddShadingSurfaces(Cell^ buildingCell, OpenStudio::Model^ osModel)
	{
		OpenStudio::ShadingSurfaceGroup^ osShadingGroup = gcnew OpenStudio::ShadingSurfaceGroup(osModel);
		List<Face^>^ faceList = buildingCell->Faces();
		int faceIndex = 0;
		for each(Face^ face in faceList)
		{
			List<Vertex^>^ vertices = face->Vertices();
			OpenStudio::Point3dVector^ facePoints = gcnew OpenStudio::Point3dVector();

			for each(Vertex^ aVertex in vertices)
			{
				Autodesk::DesignScript::Geometry::Point^ p = safe_cast<Autodesk::DesignScript::Geometry::Point^>(aVertex->Geometry);
				OpenStudio::Point3d^ aPoint = gcnew OpenStudio::Point3d(p->X, p->Y, p->Z);
				facePoints->Add(aPoint);
			}

			OpenStudio::ShadingSurface^ aShadingSurface = gcnew OpenStudio::ShadingSurface(facePoints, osModel);

			String^ surfaceName = buildingCell->ToString() + "_SHADINGSURFACE_" + (faceIndex.ToString());
			aShadingSurface->setName(surfaceName);
			aShadingSurface->setShadingSurfaceGroup(osShadingGroup);

			++faceIndex;
		}
	}

	OpenStudio::Surface^ DsosUtility::AddSurface(
		Face^ buildingFace,
		OpenStudio::Point3dVector^ osFacePoints,
		OpenStudio::Space^ osSpace,
		OpenStudio::Model^ osModel,
		Autodesk::DesignScript::Geometry::Vector^ upVector,
		double glazingRatio)
	{
		OpenStudio::Construction^ osInteriorCeilingType = nullptr;
		OpenStudio::Construction^ osExteriorRoofType = nullptr;
		OpenStudio::Construction^ osInteriorFloorType = nullptr;
		OpenStudio::Construction^ osInteriorWallType = nullptr;
		OpenStudio::Construction^ osExteriorDoorType = nullptr;
		OpenStudio::Construction^ osExteriorWallType = nullptr;
		OpenStudio::Construction^ osExteriorWindowType = nullptr;

		OpenStudio::ConstructionVector^ osConstructionTypes = osModel->getConstructions();
		OpenStudio::ConstructionVector::ConstructionVectorEnumerator^ osConstructionTypesEnumerator =
			osConstructionTypes->GetEnumerator();
		int constructionTypeCount = osConstructionTypes->Count;

		while (osConstructionTypesEnumerator->MoveNext())
		{
			OpenStudio::Construction^ osConstruction = osConstructionTypesEnumerator->Current;
			OpenStudio::OptionalString^ osConstructionTypeOptionalString = osConstruction->name();
			String^ constructionTypeName = osConstructionTypeOptionalString->__str__();
			if (constructionTypeName->Equals("000 Interior Ceiling"))
			{
				osInteriorCeilingType = osConstruction;
			}else if (constructionTypeName->Equals("000 Interior Floor"))
			{
				osInteriorFloorType = osConstruction;
			}
			else if (constructionTypeName->Equals("000 Interior Wall"))
			{
				osInteriorWallType = osConstruction;
			}
			else if (constructionTypeName->Equals("ASHRAE 189.1-2009 ExtWindow ClimateZone 4-5"))
			{
				osExteriorWindowType = osConstruction;
			}
			else if (constructionTypeName->Equals("000 Exterior Door"))
			{
				osExteriorDoorType = osConstruction;
			}
			else if (constructionTypeName->Equals("ASHRAE 189.1-2009 ExtRoof IEAD ClimateZone 2-5"))
			{
				osExteriorRoofType = osConstruction;
			}
			else if (constructionTypeName->Equals("ASHRAE 189.1-2009 ExtWall SteelFrame ClimateZone 4-8"))
			{
				osExteriorWallType = osConstruction;
			}
		} // while (osConstructionTypesEnumerator.MoveNext())

		bool isUnderground = IsUnderground(buildingFace);
		int faceType = FaceType(buildingFace, upVector);

		OpenStudio::Surface^ osSurface = gcnew OpenStudio::Surface(osFacePoints, osModel);
		osSurface->setSpace(osSpace);
		OpenStudio::OptionalString^ osSpaceOptionalString = osSpace->name();
		String^ spaceName = osSpaceOptionalString->__str__();
		String^ surfaceName = spaceName + "_SURFACE_" + buildingFace->ToString();
		osSurface->setName(surfaceName);

		int adjCount = AdjacentCellCount(buildingFace);

		if ((faceType == 1) && (adjCount > 1))
		{

			osSurface->setOutsideBoundaryCondition("Surface");
			osSurface->setSurfaceType("RoofCeiling");
			osSurface->setConstruction(osInteriorCeilingType);
			osSurface->setSunExposure("NoSun");
			osSurface->setWindExposure("NoWind");
		}
		else if ((faceType == 1) && (adjCount < 2) && (!isUnderground))
		{
			osSurface->setOutsideBoundaryCondition("Outdoors");
			osSurface->setSurfaceType("RoofCeiling");
			osSurface->setConstruction(osExteriorRoofType);
			osSurface->setSunExposure("SunExposed");
			osSurface->setWindExposure("WindExposed");
		}
		else if ((faceType == 1) && (adjCount < 2) && isUnderground)
		{
			osSurface->setOutsideBoundaryCondition("Ground");
			osSurface->setSurfaceType("RoofCeiling");
			osSurface->setConstruction(osExteriorRoofType);
			osSurface->setSunExposure("NoSun");
			osSurface->setWindExposure("NoWind");
		}
		else if ((faceType == 2) && (adjCount > 1))
		{
			osSurface->setOutsideBoundaryCondition("Surface");
			osSurface->setSurfaceType("Floor");
			osSurface->setConstruction(osInteriorFloorType);
			osSurface->setSunExposure("NoSun");
			osSurface->setWindExposure("NoWind");
		}
		else if ((faceType == 2) && (adjCount < 2))
		{
			osSurface->setOutsideBoundaryCondition("Ground");
			osSurface->setSurfaceType("Floor");
			osSurface->setConstruction(osExteriorWallType);
			osSurface->setSunExposure("NoSun");
			osSurface->setWindExposure("NoWind");
		}
		else if ((faceType == 3) && (adjCount > 1))
		{
			osSurface->setOutsideBoundaryCondition("Surface");
			osSurface->setSurfaceType("Wall");
			osSurface->setConstruction(osInteriorWallType);
			osSurface->setSunExposure("NoSun");
			osSurface->setWindExposure("NoWind");
		}
		else if ((faceType == 3) && (adjCount < 2) && isUnderground)
		{
			osSurface->setOutsideBoundaryCondition("Ground");
			osSurface->setSurfaceType("Wall");
			osSurface->setConstruction(osExteriorWallType);
			osSurface->setSunExposure("NoSun");
			osSurface->setWindExposure("NoWind");
		}
		else if ((faceType == 3) && (adjCount < 2) && (!isUnderground))
		{
			osSurface->setOutsideBoundaryCondition("Outdoors");
			osSurface->setSurfaceType("Wall");
			osSurface->setConstruction(osExteriorWallType);
			osSurface->setSunExposure("SunExposed");
			osSurface->setWindExposure("WindExposed");

			if (glazingRatio > 0)
			{
				// Triangulate the Windows
				List<Vertex^>^ scaledVertices = ScaleFaceVertices(buildingFace, glazingRatio);

				for (int i = 0; i < scaledVertices->Count - 2; ++i)
				{
					OpenStudio::Point3dVector^ osWindowFacePoints = gcnew OpenStudio::Point3dVector();

					Autodesk::DesignScript::Geometry::Point^ pDynamoPoint1 = 
						safe_cast<Autodesk::DesignScript::Geometry::Point^>(scaledVertices[0]->Geometry);
					OpenStudio::Point3d^ p1 = gcnew OpenStudio::Point3d(
						pDynamoPoint1->X,
						pDynamoPoint1->Y,
						pDynamoPoint1->Z);

					Autodesk::DesignScript::Geometry::Point^ pDynamoPoint2 =
						safe_cast<Autodesk::DesignScript::Geometry::Point^>(scaledVertices[i+1]->Geometry);
					OpenStudio::Point3d^ p2 = gcnew OpenStudio::Point3d(
						pDynamoPoint2->X,
						pDynamoPoint2->Y,
						pDynamoPoint2->Z);

					Autodesk::DesignScript::Geometry::Point^ pDynamoPoint3 =
						safe_cast<Autodesk::DesignScript::Geometry::Point^>(scaledVertices[i + 2]->Geometry);
					OpenStudio::Point3d^ p3 = gcnew OpenStudio::Point3d(
						pDynamoPoint3->X,
						pDynamoPoint3->Y,
						pDynamoPoint3->Z);

					osWindowFacePoints->Add(p1);
					osWindowFacePoints->Add(p2);
					osWindowFacePoints->Add(p3);

					OpenStudio::SubSurface^ osWindowSubSurface = gcnew OpenStudio::SubSurface(osWindowFacePoints, osModel);
					osWindowSubSurface->setSubSurfaceType("FixedWindow");
					osWindowSubSurface->setSurface(osSurface);

					double grossSubsurfaceArea = osWindowSubSurface->grossArea();
					double netSubsurfaceArea = osWindowSubSurface->netArea();
					double grossSurfaceArea = osSurface->grossArea();
					double netSurfaceArea = osSurface->netArea();
				} // for (int i = 0; i < scaledVertices->Count - 2; ++i)
			}
			else
			{
				osSurface->setWindowToWallRatio(glazingRatio, 900.0, true);
			}
		}

		return osSurface;
	}

	List<Vertex^>^ DsosUtility::ScaleFaceVertices(Face^ buildingFace, double scaleFactor)
	{
		List<Vertex^>^ scaledVertices = gcnew List<Vertex^>();
		Vertex^ faceCentre = GetFaceCentre(buildingFace);
		Autodesk::DesignScript::Geometry::Point^ faceCenterPoint =
			safe_cast<Autodesk::DesignScript::Geometry::Point^>(faceCentre->Geometry);

		List<Vertex^>^ vertices = buildingFace->Vertices();

		double sqrtScaleFactor = Math::Sqrt(scaleFactor);
		for each(Vertex^ aVertex in vertices)
		{
			Autodesk::DesignScript::Geometry::Point^ originalPoint =
				safe_cast<Autodesk::DesignScript::Geometry::Point^>(aVertex->Geometry);
			Autodesk::DesignScript::Geometry::Point^ scaledVertex = originalPoint->Subtract(faceCenterPoint->AsVector());
			scaledVertex = safe_cast<Autodesk::DesignScript::Geometry::Point^>(scaledVertex->Scale(sqrtScaleFactor));
			scaledVertex = scaledVertex->Add(faceCenterPoint->AsVector());
			scaledVertices->Add(Vertex::ByPoint(scaledVertex));
		}
		return scaledVertices;
	}

	Vertex^ DsosUtility::GetFaceCentre(Face^ buildingFace)
	{
		List<Vertex^>^ vertices = buildingFace->Vertices();
		Autodesk::DesignScript::Geometry::Point^ sumPoint = Autodesk::DesignScript::Geometry::Point::ByCoordinates(0, 0, 0);

		// assume vertices.count > 0
		if (vertices->Count < 3)
		{
			throw gcnew Exception("Invalid face");
		}

		for each(Vertex^ v in vertices)
		{
			Autodesk::DesignScript::Geometry::Point^ p =
				safe_cast<Autodesk::DesignScript::Geometry::Point^>(v->Geometry);
			sumPoint = sumPoint->Add(p->AsVector());
		}
		return Vertex::ByPoint(safe_cast<Autodesk::DesignScript::Geometry::Point^>(sumPoint->Scale(1.0 / (double)vertices->Count)));
	}

	OpenStudio::Point3dVector^ DsosUtility::GetFacePoints(Face^ buildingFace)
	{
		List<Vertex^>^ vertices = buildingFace->Vertices();
		OpenStudio::Point3dVector^ osFacePoints = gcnew OpenStudio::Point3dVector();

		for each(Vertex^ v in vertices)
		{
			Autodesk::DesignScript::Geometry::Point^ point =
				safe_cast<Autodesk::DesignScript::Geometry::Point^>(v->Geometry);
			OpenStudio::Point3d^ osPoint = gcnew OpenStudio::Point3d(
				point->X,
				point->Y,
				point->Z);

			osFacePoints->Add(osPoint);
		}

		return osFacePoints;
	}

	bool DsosUtility::IsUnderground(Face^ buildingFace)
	{
		List<Vertex^>^ vertices = buildingFace->Vertices();

		for each(Vertex^ aVertex in vertices)
		{
			Autodesk::DesignScript::Geometry::Point^ point =
				safe_cast<Autodesk::DesignScript::Geometry::Point^>(aVertex->Geometry);

			if (point->Z > 0.0)
			{
				return false;
			}
		}

		return true;
	}

	int DsosUtility::FaceType(Face^ buildingFace, Autodesk::DesignScript::Geometry::Vector^ upVector)
	{
		int returnValue = 3;
		List<Vertex^>^ vertices = buildingFace->Vertices();
		Autodesk::DesignScript::Geometry::Point^ p1 =
			safe_cast<Autodesk::DesignScript::Geometry::Point^>(vertices[0]->Geometry);
		Autodesk::DesignScript::Geometry::Point^ p2 =
			safe_cast<Autodesk::DesignScript::Geometry::Point^>(vertices[1]->Geometry);
		Autodesk::DesignScript::Geometry::Point^ p3 =
			safe_cast<Autodesk::DesignScript::Geometry::Point^>(vertices[2]->Geometry);
		Autodesk::DesignScript::Geometry::Vector^ dir = (p2->Subtract(p1->AsVector()))->AsVector()->Cross((p3->Subtract(p1->AsVector()))->AsVector());
		Autodesk::DesignScript::Geometry::Vector^ faceNormal = dir->Normalized();
		double faceAngle = faceNormal->AngleWithVector(upVector);
		if (faceAngle < 5.0)
		{
			returnValue = 1;
		}
		else if (faceAngle > 175.0)
		{
			returnValue = 2;
		}
		return returnValue;
	}

	int DsosUtility::AdjacentCellCount(Face^ buildingFace)
	{
		return buildingFace->Cells()->Count;
	}

	int DsosUtility::StoreyNumber(Cell^ buildingCell, double buildingHeight, List<double>^ floorLevels)
	{
		double volume = buildingCell->Volume();
		Vertex^ centreOfMass = buildingCell->CenterOfMass();
		Autodesk::DesignScript::Geometry::Point^ centrePoint =
			safe_cast<Autodesk::DesignScript::Geometry::Point^>(centreOfMass->Geometry);
		for (int i = 0; i < floorLevels->Count - 1; ++i)
		{
			if (centrePoint->Z > floorLevels[i] && centrePoint->Z < floorLevels[i + 1])
			{
				return i;
			}
		}

		return 0;
	}

}