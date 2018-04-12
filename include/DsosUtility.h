#pragma once

using namespace System;
using namespace System::Collections::Generic;
using namespace TopoLogic;

namespace Dsos
{
	public ref class DsosUtility
	{
	public:
		/// <summary>
		/// 
		/// </summary>
		/// <param name="osModel"></param>
		/// <param name="idfPathName"></param>
		/// <returns></returns>
		static bool CreateIdfFile(OpenStudio::Model^ osModel, String^ idfPathName);

		/// <summary>
		/// 
		/// </summary>
		/// <param name="osmTemplatePath"></param>
		/// <param name="epwWeatherPath"></param>
		/// <param name="ddyPath"></param>
		/// <param name="osmOutputPath"></param>
		/// <param name="buildingCellComplex">The building as a cell complex</param>
		/// <param name="buildingName"></param>
		/// <param name="buildingType"></param>
		/// <param name="spaceType"></param>
		/// <param name="heatingTemp"></param>
		/// <param name="coolingTemp"></param>
		/// <param name="buildingHeight"></param>
		/// <param name="floorLevels"></param>
		/// <param name="numFloors"></param>
		/// <param name="glazingRatio"></param>
		/// <param name="contextBuildings"></param>
		/// <returns></returns>
		static OpenStudio::Model^ BuildOsModel(
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
			[Autodesk::DesignScript::Runtime::DefaultArgument("null")] List<Cell^>^ contextBuildings);

		static void PerformEnergyAnalysis(String^ strOsmPath, String^ epwPathName, String^ oswPathName, String^ openStudioExePath);

		static bool SaveModel(OpenStudio::Model^ osModel, String^ osmPathName);

		static OpenStudio::Model^ GetModelFromTemplate(String^ osmTemplatePath, String^ epwWeatherPath, String^ ddyPath);

		static OpenStudio::ThermalZone^ CreateThermalZone(OpenStudio::Model^ model, OpenStudio::Space^ space, double ceilingHeight, double heatingTemp, double coolingTemp);

		static OpenStudio::BuildingStory^ AddBuildingStory(OpenStudio::Model^ model, int floorNumber);

		static OpenStudio::Building^ ComputeBuilding(
			OpenStudio::Model^ osModel,
			String^ buildingName,
			String^ buildingType,
			double buildingHeight,
			int numFloors,
			String^ spaceType);
		
		
		static OpenStudio::DefaultScheduleSet^ getDefaultScheduleSet(OpenStudio::Model^ model);

		static OpenStudio::DefaultConstructionSet^ getDefaultConstructionSet(OpenStudio::Model^ model);

		static List<OpenStudio::BuildingStory^>^ CreateBuildingStories(OpenStudio::Model^ osModel, int numFloors);

		/// <summary>
		/// 
		/// </summary>
		/// <param name="face"></param>
		/// <param name="apertureDesign"></param>
		/// <param name="numEdgeSamples"></param>
		/// <returns></returns>
		static Face^ ApplyAperture(Face^ face, Face^ apertureDesign, int numEdgeSamples);

	private:
		DsosUtility() {}

		static OpenStudio::Space^ AddSpace(
			Cell^ cell,
			OpenStudio::Model^ osModel,
			Autodesk::DesignScript::Geometry::Vector^ upVector,
			double buildingHeight,
			List<double>^ floorLevels,
			double glazingRatio,
			double heatingTemp,
			double coolingTemp
		);

		static void AddShadingSurfaces(Cell^ buildingCell, OpenStudio::Model^ osModel);
		
		static OpenStudio::Surface^ AddSurface(
			Face^ buildingFace,
			OpenStudio::Point3dVector^ osFacePoints,
			OpenStudio::Space^ osSpace,
			OpenStudio::Model^ osModel,
			Autodesk::DesignScript::Geometry::Vector^ upVector,
			double glazingRatio);

		static List<Vertex^>^ ScaleFaceVertices(Face^ buildingFace, double scaleFactor);

		static Vertex^ GetFaceCentre(Face^ buildingFace);

		static OpenStudio::Point3dVector^ GetFacePoints(Face^ buildingFace);

		static bool IsUnderground(Face^ buildingFace);

		static int FaceType(Face^ buildingFace, Autodesk::DesignScript::Geometry::Vector^ upVector);

		static int AdjacentCellCount(Face^ buildingFace);

		static int StoreyNumber(
			Cell^ buildingCell,
			double buildingHeight,
			List<double>^ floorLevels
		);

		static List<OpenStudio::BuildingStory^>^ buildingStories;
		static OpenStudio::DefaultConstructionSet^ defaultConstructionSet;
		static OpenStudio::DefaultScheduleSet^ defaultScheduleSet;

	};
}