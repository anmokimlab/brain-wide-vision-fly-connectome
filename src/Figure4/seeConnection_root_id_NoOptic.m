function [upstreamNeurons] = seeConnection_root_id_NoOptic(Want_root_ids,FAFBConnections,FAFB_consolidated_cell_types)

idx=ismember(FAFBConnections.post_root_id,Want_root_ids);
upstreamConnection=FAFBConnections(idx,:);
opticlobes={'LA_R','AME_R','ME_R','LO_R','LOP_R','LA_L','AME_L','ME_L','LO_L','LOP_L'};
TotalInSynapse=sum(upstreamConnection.syn_count);
upstreamConnection(ismember(upstreamConnection.neuropil,opticlobes),:)=[];

[upstreamNeurons,~,ic]=unique(upstreamConnection.pre_root_id);
for i=1:1:size(upstreamNeurons,1)
    idx=ic==i;
    upstreamNeurons(i,2)=sum(upstreamConnection.syn_count(idx));
end
    upstreamNeurons=sortrows(upstreamNeurons,2,'descend');

end