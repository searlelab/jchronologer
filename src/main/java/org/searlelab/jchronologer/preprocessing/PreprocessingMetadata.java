package org.searlelab.jchronologer.preprocessing;

import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.List;
import java.util.Map;

public final class PreprocessingMetadata {

    @JsonProperty("aa_to_int")
    private Map<String, Integer> aaToInt;

    @JsonProperty("mod_regex_rules")
    private List<RegexRuleSpec> modRegexRules;

    @JsonProperty("nterm_keys")
    private Map<String, String> ntermKeys;

    @JsonProperty("max_peptide_len")
    private int maxPeptideLen;

    public Map<String, Integer> getAaToInt() {
        return aaToInt;
    }

    public List<RegexRuleSpec> getModRegexRules() {
        return modRegexRules;
    }

    public Map<String, String> getNtermKeys() {
        return ntermKeys;
    }

    public int getMaxPeptideLen() {
        return maxPeptideLen;
    }

    public int tokenArrayLength() {
        return maxPeptideLen + 2;
    }

    public static final class RegexRuleSpec {
        private String pattern;
        private String token;

        public String getPattern() {
            return pattern;
        }

        public String getToken() {
            return token;
        }
    }
}
