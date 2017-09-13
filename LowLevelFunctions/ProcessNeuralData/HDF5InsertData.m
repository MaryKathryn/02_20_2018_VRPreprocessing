function HDF5InsertData(nodeID,dsetName,data)%, varargin)
%%% This is a helper function to insert data into a node in an HDF5 file
%%% tree. We first create a dataspace according to dType, then create a
%%% dataset at nodeID, then write to that dataset with our input data.

% Important! The only types accepted are: doubles, ints, strings, and

%%%Lyndon Duong Aug 2015

%%




%infer the data type if no data type
    function inferredType = InferDataType(dataToInfer)

        if ischar(dataToInfer)
            inferredType = 'H5T_C_S1';
        elseif isinteger(dataToInfer)
            inferredType='H5T_NATIVE_INT';
        elseif islogical(dataToInfer)
            inferredType='H5T_STD_U16BE';
        elseif isa(dataToInfer,'double')
            inferredType = 'H5T_NATIVE_DOUBLE'; %default to double.. beware of errors
        end
    end


inferredType = InferDataType(data); %if no datatype input, we must infer the datatype
typeID = H5T.copy(inferredType);    

% Frustratingly, Matlab cannot save logicals in HDF5 file formats (see
% https://www.hdfgroup.org/hdf5-quest.html#bool). Thus, we must convert
% Rasters to uint8 prior to saving.
% Rasters = uint8(Rasters);
if islogical(data)
    data = uint8(data);
end

chunk = fliplr(size(data));

% unlimited = H5ML.get_constant_value('H5S_UNLIMITED');
dataspaceID = H5S.create_simple(ndims(data), fliplr(size(data)),[]);%,ones(ndims(data),1)*-1);

% Create the dataset creation property list, add the gzip
% compression filter and set the chunk size.
dcpl = H5P.create('H5P_DATASET_CREATE');
H5P.set_deflate(dcpl, 4); 


H5P.set_chunk(dcpl, chunk);


datasetID = H5D.create(nodeID,dsetName,typeID,dataspaceID,dcpl);
H5D.write(datasetID,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',data);



H5P.close(dcpl);
H5T.close(typeID);
H5D.close(datasetID); %close the dataset
H5S.close(dataspaceID); %close the dataspace


end