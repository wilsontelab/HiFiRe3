//! RE site-specific behaviors during benchmarking.

/// An RE name, site and INFO tag.
pub struct ReSite {
    pub name: String,
    pub site: String,
    pub size_tag:     Vec<u8>,
    pub size_tag_str: String,
    pub dist_tag:     Vec<u8>,
    pub dist_tag_str: String,
}
impl ReSite {
    pub fn from_re_name(re_name: &str, site: Option<&str>) -> Self {
        let (size_tag, size_tag_str) = get_re_tag_name(re_name, "SIZE");
        let (dist_tag, dist_tag_str) = get_re_tag_name(re_name, "DIST");
        ReSite {
            name: re_name.to_string(),
            site: site.unwrap_or("NA").to_string(),
            size_tag,
            size_tag_str,
            dist_tag,
            dist_tag_str
        }  
    }
}
/// Get the INFO tag name for a given RE site fragment size.
fn get_re_tag_name(re_name: &str, tag_type: &str) -> (Vec<u8>, String) {
    let tag = format!("{}{}", re_name.to_ascii_uppercase(), tag_type);
    (tag.as_bytes().to_vec(), tag)
}
