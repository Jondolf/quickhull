//! This example demonstrates the built-in 3d shapes in Bevy.
//! The scene includes a patterned texture and a rotation for visualizing the normals and UVs.

use std::f32::consts::PI;

use bevy::{
    color::palettes::basic::SILVER,
    pbr::wireframe::{Wireframe, WireframePlugin},
    prelude::*,
    render::{
        mesh::{PrimitiveTopology, VertexAttributeValues},
        render_asset::RenderAssetUsages,
        render_resource::{Extent3d, TextureDimension, TextureFormat},
    },
};
use quickhull::ConvexHull;

fn main() {
    App::new()
        .add_plugins((
            DefaultPlugins.set(ImagePlugin::default_nearest()),
            WireframePlugin,
        ))
        .add_systems(Startup, setup)
        .add_systems(Update, (rotate, render_convex_hulls))
        .run();
}

/// A marker component for our shapes so we can query them separately from the ground plane
#[derive(Component)]
struct Shape;

const SHAPES_X_EXTENT: f32 = 14.0;
const EXTRUSION_X_EXTENT: f32 = 16.0;
const Z_EXTENT: f32 = 5.0;

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut images: ResMut<Assets<Image>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    let debug_material = materials.add(StandardMaterial {
        base_color_texture: Some(images.add(uv_debug_texture())),
        ..default()
    });
    let hull_material = materials.add(Color::Srgba(Srgba::new(0.7, 0.7, 0.7, 0.4)));

    let shapes = [
        meshes.add(Cuboid::default()),
        meshes.add(Tetrahedron::default()),
        meshes.add(Capsule3d::default()),
        meshes.add(Torus::default()),
        meshes.add(Cylinder::default()),
        meshes.add(Cone::default()),
        meshes.add(ConicalFrustum::default()),
        meshes.add(Sphere::default().mesh().ico(5).unwrap()),
        meshes.add(Sphere::default().mesh().uv(32, 18)),
    ];

    let extrusions = [
        meshes.add(Extrusion::new(Rectangle::default(), 1.)),
        meshes.add(Extrusion::new(Capsule2d::default(), 1.)),
        meshes.add(Extrusion::new(Annulus::default(), 1.)),
        meshes.add(Extrusion::new(Circle::default(), 1.)),
        meshes.add(Extrusion::new(Ellipse::default(), 1.)),
        meshes.add(Extrusion::new(RegularPolygon::default(), 1.)),
        meshes.add(Extrusion::new(Triangle2d::default(), 1.)),
    ];

    let num_shapes = shapes.len();

    for (i, shape) in shapes.into_iter().enumerate() {
        commands
            .spawn((
                PbrBundle {
                    mesh: shape,
                    material: debug_material.clone(),
                    transform: Transform::from_xyz(
                        -SHAPES_X_EXTENT / 2.
                            + i as f32 / (num_shapes - 1) as f32 * SHAPES_X_EXTENT,
                        2.0,
                        Z_EXTENT / 2.,
                    )
                    .with_rotation(Quat::from_rotation_x(-PI / 4.)),
                    ..default()
                },
                Shape,
            ))
            .with_children(|commands| {
                commands.spawn((
                    PbrBundle {
                        mesh: meshes.add(Mesh::new(
                            PrimitiveTopology::default(),
                            RenderAssetUsages::default(),
                        )),
                        transform: Transform::from_scale(Vec3::splat(1.1)),
                        material: hull_material.clone(),
                        ..default()
                    },
                    Wireframe,
                ));
            });
    }

    let num_extrusions = extrusions.len();

    for (i, shape) in extrusions.into_iter().enumerate() {
        commands
            .spawn((
                PbrBundle {
                    mesh: shape,
                    material: debug_material.clone(),
                    transform: Transform::from_xyz(
                        -EXTRUSION_X_EXTENT / 2.
                            + i as f32 / (num_extrusions - 1) as f32 * EXTRUSION_X_EXTENT,
                        2.0,
                        -Z_EXTENT / 2.,
                    )
                    .with_rotation(Quat::from_rotation_x(-PI / 4.)),
                    ..default()
                },
                Shape,
            ))
            .with_children(|commands| {
                commands.spawn((
                    PbrBundle {
                        mesh: meshes.add(Mesh::new(
                            PrimitiveTopology::default(),
                            RenderAssetUsages::default(),
                        )),
                        material: hull_material.clone(),
                        transform: Transform::from_scale(Vec3::splat(1.1)),
                        ..default()
                    },
                    Wireframe,
                ));
            });
    }

    commands.spawn(PointLightBundle {
        point_light: PointLight {
            shadows_enabled: true,
            intensity: 10_000_000.,
            range: 100.0,
            shadow_depth_bias: 0.2,
            ..default()
        },
        transform: Transform::from_xyz(8.0, 16.0, 8.0),
        ..default()
    });

    // ground plane
    commands.spawn(PbrBundle {
        mesh: meshes.add(Plane3d::default().mesh().size(50.0, 50.0).subdivisions(10)),
        material: materials.add(Color::from(SILVER)),
        ..default()
    });

    commands.spawn(Camera3dBundle {
        transform: Transform::from_xyz(0.0, 7., 14.0).looking_at(Vec3::new(0., 1., 0.), Vec3::Y),
        ..default()
    });

    /*commands.spawn(
        TextBundle::from_section("Press space to toggle wireframes", TextStyle::default())
            .with_style(Style {
                position_type: PositionType::Absolute,
                top: Val::Px(12.0),
                left: Val::Px(12.0),
                ..default()
            }),
    );*/
}

fn rotate(mut query: Query<&mut Transform, With<Shape>>, time: Res<Time>) {
    for mut transform in &mut query {
        transform.rotate_y(time.delta_seconds() / 2.);
    }
}

/// Creates a colorful test pattern
fn uv_debug_texture() -> Image {
    const TEXTURE_SIZE: usize = 8;

    let mut palette: [u8; 32] = [
        255, 102, 159, 255, 255, 159, 102, 255, 236, 255, 102, 255, 121, 255, 102, 255, 102, 255,
        198, 255, 102, 198, 255, 255, 121, 102, 255, 255, 236, 102, 255, 255,
    ];

    let mut texture_data = [0; TEXTURE_SIZE * TEXTURE_SIZE * 4];
    for y in 0..TEXTURE_SIZE {
        let offset = TEXTURE_SIZE * y * 4;
        texture_data[offset..(offset + TEXTURE_SIZE * 4)].copy_from_slice(&palette);
        palette.rotate_right(4);
    }

    Image::new_fill(
        Extent3d {
            width: TEXTURE_SIZE as u32,
            height: TEXTURE_SIZE as u32,
            depth_or_array_layers: 1,
        },
        TextureDimension::D2,
        &texture_data,
        TextureFormat::Rgba8UnormSrgb,
        RenderAssetUsages::RENDER_WORLD,
    )
}

fn render_convex_hulls(
    query: Query<(&Handle<Mesh>, &Children)>,
    hull_query: Query<&Handle<Mesh>>,
    mut meshes: ResMut<Assets<Mesh>>,
) {
    for (handle, children) in &query {
        let Some(mesh) = meshes.get(handle) else {
            continue;
        };

        let Some(VertexAttributeValues::Float32x3(positions)) =
            mesh.attribute(Mesh::ATTRIBUTE_POSITION)
        else {
            continue;
        };

        let positions = positions
            .iter()
            .map(|v| Vec3::from(*v).as_dvec3())
            .collect::<Vec<_>>();

        let Ok(hull) = ConvexHull::try_new(&positions, None).map_err(|e| warn!("{:?}", e)) else {
            continue;
        };

        for handle in hull_query.iter_many(children) {
            let Some(mesh) = meshes.get_mut(handle) else {
                continue;
            };

            let (vertices, indices) = hull.vertices_indices();

            mesh.insert_attribute(
                Mesh::ATTRIBUTE_POSITION,
                vertices
                    .iter()
                    .map(|v| v.as_vec3().to_array())
                    .collect::<Vec<_>>(),
            );

            mesh.insert_indices(bevy::render::mesh::Indices::U32(
                indices.iter().map(|i| *i as u32).collect::<Vec<_>>(),
            ));

            mesh.duplicate_vertices();

            mesh.compute_flat_normals();
        }
    }
}
