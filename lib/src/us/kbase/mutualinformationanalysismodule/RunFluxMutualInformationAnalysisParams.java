
package us.kbase.mutualinformationanalysismodule;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: RunFluxMutualInformationAnalysisParams</p>
 * 
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "fbamodel_id",
    "compounds",
    "workspace",
    "media_id"
})
public class RunFluxMutualInformationAnalysisParams {

    @JsonProperty("fbamodel_id")
    private java.lang.String fbamodelId;
    @JsonProperty("compounds")
    private List<String> compounds;
    @JsonProperty("workspace")
    private java.lang.String workspace;
    @JsonProperty("media_id")
    private java.lang.String mediaId;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("fbamodel_id")
    public java.lang.String getFbamodelId() {
        return fbamodelId;
    }

    @JsonProperty("fbamodel_id")
    public void setFbamodelId(java.lang.String fbamodelId) {
        this.fbamodelId = fbamodelId;
    }

    public RunFluxMutualInformationAnalysisParams withFbamodelId(java.lang.String fbamodelId) {
        this.fbamodelId = fbamodelId;
        return this;
    }

    @JsonProperty("compounds")
    public List<String> getCompounds() {
        return compounds;
    }

    @JsonProperty("compounds")
    public void setCompounds(List<String> compounds) {
        this.compounds = compounds;
    }

    public RunFluxMutualInformationAnalysisParams withCompounds(List<String> compounds) {
        this.compounds = compounds;
        return this;
    }

    @JsonProperty("workspace")
    public java.lang.String getWorkspace() {
        return workspace;
    }

    @JsonProperty("workspace")
    public void setWorkspace(java.lang.String workspace) {
        this.workspace = workspace;
    }

    public RunFluxMutualInformationAnalysisParams withWorkspace(java.lang.String workspace) {
        this.workspace = workspace;
        return this;
    }

    @JsonProperty("media_id")
    public java.lang.String getMediaId() {
        return mediaId;
    }

    @JsonProperty("media_id")
    public void setMediaId(java.lang.String mediaId) {
        this.mediaId = mediaId;
    }

    public RunFluxMutualInformationAnalysisParams withMediaId(java.lang.String mediaId) {
        this.mediaId = mediaId;
        return this;
    }

    @JsonAnyGetter
    public Map<java.lang.String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(java.lang.String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public java.lang.String toString() {
        return ((((((((((("RunFluxMutualInformationAnalysisParams"+" [fbamodelId=")+ fbamodelId)+", compounds=")+ compounds)+", workspace=")+ workspace)+", mediaId=")+ mediaId)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
