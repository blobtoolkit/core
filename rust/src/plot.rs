//!
//! Invoked by calling:
//! `blobtk plot <args>`

use std::collections::HashMap;

use crate::blobdir;
use crate::cli;
use crate::plot::blob::BlobData;
// use crate::io;

pub use cli::PlotOptions;
use colorous;
use svg::Document;
use usvg::{fontdb, TreeParsing, TreeTextToPath};

use self::blob::BlobDimensions;

/// Blob plot functions.
pub mod blob;

/// Category functions.
pub mod category;

/// Chart components.
pub mod component;

/// Scatter plot functions.
pub mod scatter;

/// Snail plot functions.
pub mod snail;

/// SVG styling functions.
pub mod style;

pub fn save_svg(document: &Document, options: &PlotOptions) {
    svg::save(options.output.as_str(), document).unwrap();
}

pub fn save_png(document: &Document, _: &PlotOptions) {
    let mut fontdb = fontdb::Database::new();
    fontdb.load_system_fonts();
    let mut buf = Vec::new();
    svg::write(&mut buf, document).unwrap();
    let opt = usvg::Options::default();
    let mut tree = usvg::Tree::from_data(&buf.as_slice(), &opt).unwrap();
    tree.convert_text(&fontdb);

    let pixmap_size = tree.size.to_screen_size();
    let mut pixmap = tiny_skia::Pixmap::new(pixmap_size.width(), pixmap_size.height()).unwrap();
    resvg::render(
        &tree,
        resvg::FitTo::Original,
        tiny_skia::Transform::default(),
        pixmap.as_mut(),
    )
    .unwrap();
    pixmap.save_png("test.png").unwrap();
}

pub fn plot_snail(meta: &blobdir::Meta, options: &cli::PlotOptions) {
    // let busco_list = meta.busco_list.clone().unwrap();
    let busco_field = meta.busco_list.clone().unwrap()[0].clone();
    let busco_values = blobdir::parse_field_busco(busco_field.0, &options.blobdir).unwrap();
    let busco_total = busco_field.1;
    let busco_lineage = busco_field.2;
    let gc_values = blobdir::parse_field_float("gc".to_string(), &options.blobdir).unwrap();
    let length_values = blobdir::parse_field_int("length".to_string(), &options.blobdir).unwrap();
    let n_values = blobdir::parse_field_float("n".to_string(), &options.blobdir);
    let ncount_values = blobdir::parse_field_int("ncount".to_string(), &options.blobdir).unwrap();
    let id = meta.id.clone();
    let record_type = meta.record_type.clone();

    let filters = blobdir::parse_filters(&options.filter);
    let wanted_indices = blobdir::set_filters(filters, &meta, &options.blobdir);

    let gc_filtered = blobdir::apply_filter_float(&gc_values, &wanted_indices);
    let n_filtered = match n_values {
        None => None,
        Some(values) => Some(blobdir::apply_filter_float(&values, &wanted_indices)),
    };
    let length_filtered = blobdir::apply_filter_int(&length_values, &wanted_indices);
    let ncount_filtered = blobdir::apply_filter_int(&ncount_values, &wanted_indices);
    let busco_filtered = blobdir::apply_filter_busco(&busco_values, &wanted_indices);

    let snail_stats = snail::snail_stats(
        &length_filtered,
        &gc_filtered,
        &n_filtered,
        &ncount_filtered,
        &busco_filtered,
        busco_total,
        busco_lineage,
        id,
        record_type,
        &options,
    );
    let document: Document = snail::svg(&snail_stats, &options);

    save_svg(&document, &options);

    save_png(&document, &options);

    // let cats = blobdir::parse_field_cat("buscogenes_family".to_string(), &options);
    // let identifiers = blobdir::parse_field_string("identifiers".to_string(), &options);
}

pub fn color_to_hex(color: colorous::Color) -> String {
    format!("#{:x}{:x}{:x}", color.r, color.g, color.b)
}

pub fn default_palette(count: usize) -> Vec<String> {
    let gradient = colorous::PAIRED;
    let mut list = vec![];
    for i in 0..count {
        let mut j = if i % 2 == 1 { i - 1 } else { i + 1 };
        j = j % 12;
        list.push(color_to_hex(gradient[j]));
    }
    list
}

pub fn set_palette(
    name: &Option<cli::Palette>,
    colors: &Option<Vec<String>>,
    count: usize,
) -> Vec<String> {
    let mut color_list = match name {
        Some(cli::Palette::Default) => default_palette(count),
        Some(cli::Palette::Paired) => {
            let gradient = colorous::PAIRED;
            let mut list = vec![];
            for i in 0..count {
                let j = i % 12;
                list.push(color_to_hex(gradient[j]));
            }
            list
        }
        Some(cli::Palette::Viridis) => {
            let gradient = colorous::VIRIDIS;
            (0..count)
                .map(|i| color_to_hex(gradient.eval_rational(i, count)))
                .collect()
        }
        None => default_palette(count),
    };
    if colors.is_some() {
        for color in colors.clone().unwrap() {
            let (index, hex) = color.split_once("=").unwrap();
            let i: usize = index.parse().unwrap();
            let mut hexcode = hex.to_string();
            if i <= count {
                hexcode = hexcode.replace("hex", "#");
                if !hexcode.starts_with("#") {
                    hexcode = format!("#{}", hexcode);
                }
                color_list[i] = hexcode;
            }
        }
    }
    color_list
}

pub fn plot_blob(meta: &blobdir::Meta, options: &cli::PlotOptions) {
    // let busco_list = meta.busco_list.clone().unwrap();
    let mut plot_meta: HashMap<String, String> = HashMap::new();
    if options.x_field.is_some() {
        plot_meta.insert("x".to_string(), options.x_field.clone().unwrap());
    } else {
        plot_meta.insert("x".to_string(), meta.plot.x.clone().unwrap());
    }
    if options.y_field.is_some() {
        plot_meta.insert("y".to_string(), options.y_field.clone().unwrap());
    } else {
        plot_meta.insert("y".to_string(), meta.plot.y.clone().unwrap());
    }
    if options.z_field.is_some() {
        plot_meta.insert("z".to_string(), options.z_field.clone().unwrap());
    } else {
        plot_meta.insert("z".to_string(), meta.plot.z.clone().unwrap());
    }
    if options.cat_field.is_some() {
        plot_meta.insert("cat".to_string(), options.cat_field.clone().unwrap());
    } else {
        plot_meta.insert("cat".to_string(), meta.plot.cat.clone().unwrap());
    }
    // TODO: handle empty values

    let (plot_values, cat_values) = blobdir::get_plot_values(&meta, &options.blobdir, &plot_meta);

    let palette = set_palette(&options.palette, &options.color, options.cat_count);

    let (cat_order, cat_indices) = category::set_cat_order(
        &cat_values,
        &options.cat_order,
        &options.cat_count,
        &palette,
    );
    // let id = meta.id.clone();
    // let record_type = meta.record_type.clone();

    let filters = blobdir::parse_filters(&options.filter);
    let wanted_indices = blobdir::set_filters(filters, &meta, &options.blobdir);
    let blob_data = BlobData {
        x: blobdir::apply_filter_float(&plot_values["x"], &wanted_indices),
        y: blobdir::apply_filter_float(&plot_values["y"], &wanted_indices),
        z: blobdir::apply_filter_float(&plot_values["z"], &wanted_indices),
        cat: blobdir::apply_filter_int(&cat_indices, &wanted_indices),
        cat_order,
    };

    let scatter_data = blob::blob_points(plot_meta, blob_data, &meta, &options);

    let dimensions = BlobDimensions {
        ..Default::default()
    };

    let document: Document = blob::svg(&dimensions, &scatter_data, &options);

    save_svg(&document, &options);

    save_png(&document, &options);
}

/// Execute the `plot` subcommand from `blobtk`.
pub fn plot(options: &cli::PlotOptions) -> Result<(), Box<dyn std::error::Error>> {
    let meta = blobdir::parse_blobdir(&options.blobdir);
    let view = &options.view;
    match view {
        Some(cli::View::Blob) => plot_blob(&meta, &options),
        Some(cli::View::Snail) => plot_snail(&meta, &options),
        _ => (),
    }
    Ok(())
}
