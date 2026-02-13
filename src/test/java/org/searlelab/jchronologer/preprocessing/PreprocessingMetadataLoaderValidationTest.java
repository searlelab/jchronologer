package org.searlelab.jchronologer.preprocessing;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.fasterxml.jackson.databind.ObjectMapper;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.preprocessing.PreprocessingMetadataLoader.CompiledPreprocessingMetadata;

class PreprocessingMetadataLoaderValidationTest {

    @Test
    void loadFromClasspathThrowsForMissingResource() {
        IllegalArgumentException error = assertThrows(
                IllegalArgumentException.class,
                () -> PreprocessingMetadataLoader.loadFromClasspath("missing/preprocessing.json"));
        assertTrue(error.getMessage().contains("Missing preprocessing resource"));
    }

    @Test
    void loadFromClasspathValidatesAaToInt() {
        IllegalArgumentException error = assertThrows(
                IllegalArgumentException.class,
                () -> PreprocessingMetadataLoader.loadFromClasspath("data/preprocessing/invalid_aa_to_int.json"));
        assertTrue(error.getMessage().contains("Invalid aa_to_int"));
    }

    @Test
    void loadFromClasspathValidatesModRegexRules() {
        IllegalArgumentException error = assertThrows(
                IllegalArgumentException.class,
                () -> PreprocessingMetadataLoader.loadFromClasspath("data/preprocessing/invalid_mod_regex_rules.json"));
        assertTrue(error.getMessage().contains("Missing mod_regex_rules"));
    }

    @Test
    void loadFromClasspathValidatesNullModRegexRuleEntry() {
        IllegalArgumentException error = assertThrows(
                IllegalArgumentException.class,
                () -> PreprocessingMetadataLoader.loadFromClasspath(
                        "data/preprocessing/invalid_mod_regex_rule_entry_null.json"));
        assertTrue(error.getMessage().contains("Invalid mod_regex_rules entry at index 0"));
    }

    @Test
    void loadFromClasspathValidatesBlankModRegexPattern() {
        IllegalArgumentException error = assertThrows(
                IllegalArgumentException.class,
                () -> PreprocessingMetadataLoader.loadFromClasspath(
                        "data/preprocessing/invalid_mod_regex_pattern_blank.json"));
        assertTrue(error.getMessage().contains("Invalid mod_regex_rules pattern at index 0"));
    }

    @Test
    void loadFromClasspathValidatesBlankModRegexToken() {
        IllegalArgumentException error = assertThrows(
                IllegalArgumentException.class,
                () -> PreprocessingMetadataLoader.loadFromClasspath(
                        "data/preprocessing/invalid_mod_regex_token_blank.json"));
        assertTrue(error.getMessage().contains("Invalid mod_regex_rules token at index 0"));
    }

    @Test
    void loadFromClasspathValidatesRegexPatternSyntax() {
        IllegalArgumentException error = assertThrows(
                IllegalArgumentException.class,
                () -> PreprocessingMetadataLoader.loadFromClasspath(
                        "data/preprocessing/invalid_mod_regex_pattern_syntax.json"));
        assertTrue(error.getMessage().contains("Invalid mod_regex_rules regex at index 0"));
    }

    @Test
    void loadFromClasspathValidatesNtermKeys() {
        IllegalArgumentException error = assertThrows(
                IllegalArgumentException.class,
                () -> PreprocessingMetadataLoader.loadFromClasspath("data/preprocessing/invalid_nterm_keys.json"));
        assertTrue(error.getMessage().contains("Missing nterm_keys"));
    }

    @Test
    void loadFromClasspathValidatesMaxPeptideLength() {
        IllegalArgumentException error = assertThrows(
                IllegalArgumentException.class,
                () -> PreprocessingMetadataLoader.loadFromClasspath("data/preprocessing/invalid_max_peptide_len.json"));
        assertTrue(error.getMessage().contains("Invalid max_peptide_len"));
    }

    @Test
    void loadFromClasspathCompilesMinimalMetadata() {
        CompiledPreprocessingMetadata metadata =
                PreprocessingMetadataLoader.loadFromClasspath("data/preprocessing/minimal_valid.json");
        assertEquals(1, metadata.getAaToInt().size());
        assertEquals(1, metadata.getModRegexRules().size());
        assertEquals("[+42.010565]", metadata.getNtermKeys().get("[+42.010565]"));
        assertEquals(6, metadata.tokenArrayLength());
    }

    @Test
    void preprocessingMetadataTokenArrayLengthMatchesMaxPlusTwo() throws Exception {
        String json = "{\"aa_to_int\":{\"A\":1},\"mod_regex_rules\":[],\"nterm_keys\":{\"N\":\"N\"},\"max_peptide_len\":9}";
        PreprocessingMetadata metadata = new ObjectMapper().readValue(json, PreprocessingMetadata.class);
        assertEquals(11, metadata.tokenArrayLength());
    }
}
