//! Type aliases for expressive code that is clear about data types in use.

/// The bit representation of a single bit-encoded sample.
pub type SampleBit = u32;
/// The bit representation of a multiple bit-encoded samples.
pub type SampleBits = u32;

/// The String name of a reference chromosome.
pub type ChromName  = String;
/// The 1-based index of a reference chromosome.
pub type ChromIndex1 = u8;
/// An integer capable of holding values up to the length of any chromosome.
pub type ChromLength = u32;

/// A 0-based chromosome position, e.g., as used for half-open span starts.
pub type ChromPos0 = u32;
/// A 1-based chromosome position, e.g., as used for half-open (i.e., 0-based
/// exclusive) span ends.
pub type ChromPos1 = u32;
/// A 0-based position on any type of DNA sequence, e.g., as used for half-open 
/// span starts.
pub type SeqPos0   = u32;
// /// A 1-based position on any type of DNA sequence, e.g., as used for half-open 
// /// (i.e., 0-based exclusive) span ends.
// pub type SeqPos1   = u32;

/// A String expected to contain only valid UTF-8 single-byte characters that
/// represents the name of a query sequence in a BAM or similar sequence file.
pub type QName = String;

/// The 0-based index of a reads among a set of related sequencing reads.
pub type ReadIndex = usize;

/// A String represention of a base sequence expected to only containing bases
/// A, C, G, T, and N in uppercase.
pub type UppercaseACGTN = String;
// /// A String represention of a base sequence expected to only containing bases
// /// A, C, G, T, and N in lowercase.
// pub type LowercaseACGTN = String;
// /// A String represention of a base sequence expected to only containing bases
// /// A, C, G, T, and N in any combination of uppercase and lowercase.
// pub type MixedcaseACGTN = String;
/// A single valid UTF-8 character byte expected to be either A, C, G, T, or N.
pub type BaseByteACGTN = u8;
/// A single valid UTF-8 character byte expected to be either A, C, G, or T only.
pub type BaseByteACGT = u8;

/// A PHRED-scaled probability that a called base is wrong without any offsets, 
/// hence the value can be used directly without subtracting 33.
pub type PhredQual = u8;

/// A comma-delimited list of values as a String, regardless of the original
/// data type.
pub type CommaDelimited = String;

// /// A boolean cast as a u8 integer for printing to file.
// pub type IntegerBool = u8;
