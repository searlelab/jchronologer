package org.searlelab.jchronologer.preprocessing;

import com.fasterxml.jackson.databind.DeserializationFeature;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

public final class PreprocessingMetadataLoader {

    private static final ObjectMapper MAPPER = new ObjectMapper()
            .configure(DeserializationFeature.FAIL_ON_UNKNOWN_PROPERTIES, false);

    private PreprocessingMetadataLoader() {
    }

    public static CompiledPreprocessingMetadata loadFromClasspath(String resource) {
        try (InputStream stream = Thread.currentThread().getContextClassLoader().getResourceAsStream(resource)) {
            if (stream == null) {
                throw new IllegalArgumentException("Missing preprocessing resource: " + resource);
            }
            PreprocessingMetadata metadata = MAPPER.readValue(stream, PreprocessingMetadata.class);
            validate(metadata, resource);
            return compile(metadata, resource);
        } catch (IOException e) {
            throw new IllegalStateException("Failed to load preprocessing resource: " + resource, e);
        }
    }

    private static void validate(PreprocessingMetadata metadata, String resource) {
        if (metadata.getAaToInt() == null || metadata.getAaToInt().isEmpty()) {
            throw new IllegalArgumentException("Invalid aa_to_int in preprocessing resource: " + resource);
        }
        if (metadata.getModRegexRules() == null) {
            throw new IllegalArgumentException("Missing mod_regex_rules in preprocessing resource: " + resource);
        }
        for (int i = 0; i < metadata.getModRegexRules().size(); i++) {
            PreprocessingMetadata.RegexRuleSpec rule = metadata.getModRegexRules().get(i);
            if (rule == null) {
                throw new IllegalArgumentException(
                        "Invalid mod_regex_rules entry at index " + i + " in preprocessing resource: " + resource);
            }
            if (rule.getPattern() == null || rule.getPattern().isBlank()) {
                throw new IllegalArgumentException(
                        "Invalid mod_regex_rules pattern at index " + i + " in preprocessing resource: " + resource);
            }
            if (rule.getToken() == null || rule.getToken().isBlank()) {
                throw new IllegalArgumentException(
                        "Invalid mod_regex_rules token at index " + i + " in preprocessing resource: " + resource);
            }
        }
        if (metadata.getNtermKeys() == null || metadata.getNtermKeys().isEmpty()) {
            throw new IllegalArgumentException("Missing nterm_keys in preprocessing resource: " + resource);
        }
        if (metadata.getMaxPeptideLen() <= 0) {
            throw new IllegalArgumentException("Invalid max_peptide_len in preprocessing resource: " + resource);
        }
    }

    private static CompiledPreprocessingMetadata compile(PreprocessingMetadata metadata, String resource) {
        List<CompiledPreprocessingMetadata.CompiledRegexRule> regexRules = new ArrayList<>();
        for (int i = 0; i < metadata.getModRegexRules().size(); i++) {
            PreprocessingMetadata.RegexRuleSpec rule = metadata.getModRegexRules().get(i);
            String normalizedPattern = normalizePythonQuantifiers(rule.getPattern());
            try {
                regexRules.add(new CompiledPreprocessingMetadata.CompiledRegexRule(
                        Pattern.compile(normalizedPattern),
                        rule.getToken()));
            } catch (PatternSyntaxException e) {
                throw new IllegalArgumentException(
                        "Invalid mod_regex_rules regex at index " + i + " in preprocessing resource: " + resource,
                        e);
            }
        }
        return new CompiledPreprocessingMetadata(
                metadata.getAaToInt(),
                regexRules,
                metadata.getNtermKeys(),
                metadata.getMaxPeptideLen());
    }

    static String normalizePythonQuantifiers(String pattern) {
        return pattern.replaceAll("\\{,(\\d+)\\}", "{0,$1}");
    }

    public static final class CompiledPreprocessingMetadata {

        private final Map<String, Integer> aaToInt;
        private final List<CompiledRegexRule> modRegexRules;
        private final Map<String, String> ntermKeys;
        private final int maxPeptideLen;

        public CompiledPreprocessingMetadata(
                Map<String, Integer> aaToInt,
                List<CompiledRegexRule> modRegexRules,
                Map<String, String> ntermKeys,
                int maxPeptideLen) {
            this.aaToInt = aaToInt;
            this.modRegexRules = modRegexRules;
            this.ntermKeys = ntermKeys;
            this.maxPeptideLen = maxPeptideLen;
        }

        public Map<String, Integer> getAaToInt() {
            return aaToInt;
        }

        public List<CompiledRegexRule> getModRegexRules() {
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

        public static final class CompiledRegexRule {
            private final Pattern pattern;
            private final String token;

            public CompiledRegexRule(Pattern pattern, String token) {
                this.pattern = pattern;
                this.token = token;
            }

            public Pattern getPattern() {
                return pattern;
            }

            public String getToken() {
                return token;
            }
        }
    }
}
