//! Demonstrates the coplanar face merging feature of [`ConvexHull3d`] for the Utah teapot model.

use bevy::{
    asset::RenderAssetUsages, camera::Exposure, core_pipeline::tonemapping::Tonemapping,
    mesh::PrimitiveTopology, prelude::*,
};
use quickhull::{ConvexHull3d, ConvexHull3dSettings};

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .init_resource::<HullSettings>()
        .add_systems(Startup, setup)
        .add_systems(
            Update,
            (change_coplanarity_tolerance, rotate, render_hull).chain(),
        )
        .run();
}

#[derive(Component)]
struct Rotating;

#[derive(Component, Default, Deref, DerefMut)]
struct PointCloud(Vec<Vec3A>);

#[derive(Component)]
struct Hull;

#[derive(Resource, Default, Deref, DerefMut)]
struct HullSettings(ConvexHull3dSettings);

#[derive(Component)]
struct CoplanarityToleranceText;

fn setup(mut commands: Commands, mut meshes: ResMut<Assets<Mesh>>, assets: Res<AssetServer>) {
    // Spawn the Utah teapot
    commands.spawn((
        SceneRoot(assets.load(GltfAssetLabel::Scene(0).from_asset("utah_teapot.glb"))),
        Transform::from_xyz(-2.5, 0.0, 0.0).with_scale(Vec3::splat(1.0)),
        Rotating,
        // This will be filled by the `on_scene_ready` observer
        PointCloud::default(),
    ));

    // Spawn the hull mesh (initially empty, will be updated by a system)
    commands.spawn((
        Mesh3d(meshes.add(Mesh::new(
            PrimitiveTopology::TriangleList,
            RenderAssetUsages::default(),
        ))),
        Transform::from_xyz(2.5, 0.0, 0.0),
        Rotating,
        Hull,
    ));

    // UI
    commands.spawn((
        CoplanarityToleranceText,
        Text::new("Coplanarity dot tolerance: 1.0"),
        Node {
            position_type: PositionType::Absolute,
            top: Val::Px(10.0),
            left: Val::Px(10.0),
            ..default()
        },
    ));
    commands.spawn((
        Text::new("Use Up/Down arrows to adjust coplanarity tolerance"),
        Node {
            position_type: PositionType::Absolute,
            top: Val::Px(40.0),
            left: Val::Px(10.0),
            ..default()
        },
    ));

    // Light
    commands.spawn((
        PointLight {
            intensity: 50_000_000.0,
            range: 100.0,
            shadows_enabled: true,
            ..default()
        },
        Transform::from_xyz(1.0, 6.0, 7.0),
    ));

    // Camera
    commands.spawn((
        Camera3d::default(),
        AmbientLight {
            brightness: 4000.0,
            color: Color::WHITE,
            ..default()
        },
        Exposure::SUNLIGHT,
        Tonemapping::AcesFitted,
        Transform::from_xyz(0.0, 2.0, 8.0).looking_at(Vec3::new(0.0, 0.75, 0.0), Vec3::Y),
    ));
}

fn rotate(mut query: Query<&mut Transform, With<Rotating>>, time: Res<Time>) {
    let angle = time.elapsed_secs() * core::f32::consts::PI / 4.0;
    for mut transform in &mut query {
        transform.rotation = Quat::from_rotation_y(angle);
    }
}

fn render_hull(
    point_cloud: Single<(Entity, &mut PointCloud)>,
    hull_transform: Single<&GlobalTransform, With<Hull>>,
    child_query: Query<&Children>,
    mesh_query: Query<&Mesh3d>,
    settings: Res<HullSettings>,
    meshes: ResMut<Assets<Mesh>>,
    mut gizmos: Gizmos,
) {
    let (point_cloud_entity, mut point_cloud) = point_cloud.into_inner();

    // Get all vertices from the loaded scene.
    point_cloud.clear();
    for mesh in mesh_query.iter_many(child_query.iter_descendants(point_cloud_entity)) {
        let mesh = meshes.get(mesh.id()).unwrap();
        let vertex_positions = mesh.attribute(Mesh::ATTRIBUTE_POSITION).unwrap();
        let vertex_positions = vertex_positions.as_float3().unwrap();
        point_cloud.extend(
            vertex_positions
                .iter()
                .map(|&[x, y, z]| Vec3A::new(x, y, z)),
        );
    }

    // Compute the convex hull.
    let Ok(hull) = ConvexHull3d::try_from_points(&point_cloud.0, settings.0) else {
        error!("Failed to compute convex hull");
        return;
    };
    let (vertices, indices, faces, planes) = hull.into_parts();

    // Render the edges of the hull using gizmos.
    for (i, face) in faces.iter().enumerate() {
        let indices = face.vertex_indices(&indices);

        // Slightly offset the edges along the face normal
        // to display them on top of the mesh.
        let offset = planes[i].normal() * 0.01;
        for j in 0..indices.len() {
            let a = vertices[indices[j] as usize] + offset;
            let b = vertices[indices[(j + 1) % indices.len()] as usize] + offset;
            gizmos.line(
                hull_transform.transform_point(Vec3::from(a.to_array())),
                hull_transform.transform_point(Vec3::from(b.to_array())),
                // Color based on index to differentiate faces
                Color::hsl(i as f32 / faces.len() as f32 * 360.0, 0.7, 0.5),
            );
        }
    }
}

fn change_coplanarity_tolerance(
    mut settings: ResMut<HullSettings>,
    mut text: Single<&mut Text, With<CoplanarityToleranceText>>,
    input: Res<ButtonInput<KeyCode>>,
    time: Res<Time>,
) {
    let change_rate = 0.2;
    let delta = change_rate * time.delta_secs();
    let acos_tolerance = settings.coplanarity_dot_tolerance.acos();
    if input.pressed(KeyCode::ArrowUp) {
        settings.coplanarity_dot_tolerance = (acos_tolerance - delta).cos().min(1.0);
    }
    if input.pressed(KeyCode::ArrowDown) {
        settings.coplanarity_dot_tolerance = (acos_tolerance + delta).cos().max(0.0);
    }
    text.0 = format!(
        "Coplanarity dot tolerance: {:.3}",
        settings.coplanarity_dot_tolerance
    );
}
